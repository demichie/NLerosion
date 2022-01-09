def advance_time( h_init , h_old , delta_t , delta_x , delta_y , lambda_wb , k_wb ,
                  S_c , enne , k , max_nlc , max_inner_iter , res ,
                  Dirichlet, Transient , Neumann , grow_rate , A_c ,
                  simtime , verbose_level ):

    # ADVANCE_TIME This function update of one time step delta_t the solution
    # of dh/dt = div( f(h) grad(h) )
    #
    # where
    #
    # f(h) = k * ( 1 + ( max_nlc -1 ) * ( |grad(h)| / S_c )^n )
    #
    # k is the diffusivity constant, S_c is the tangent of the critical slope.
    # and max_nlc is the maximum value of the non-linear coefficient that
    # multiply k.

    import numpy as np
    from nonlinear_function import nonlinear_function
    from numerical_fluxes import numerical_fluxes
    from relative_coeff import relative_coeff

    nx = h_init.shape[0]
    ny = h_init.shape[1]

    h_temp = np.zeros((nx,ny))

    h_temp[0:nx,0:ny] = h_init[0:nx,0:ny]
    h_new = np.zeros((nx,ny))

    w = np.zeros((nx,ny)) # temporary solution array
    v = np.zeros((nx,ny)) # temporary solution array

    ax = np.zeros((nx)) # lower-diagonal
    bx = np.zeros((nx)) # diagonal
    cx = np.zeros((nx)) # upper-diagonal
    dx = np.zeros((nx)) # right-hand side

    ay = np.zeros((ny)) # lower-diagonal
    by = np.zeros((ny)) # diagonal
    cy = np.zeros((ny)) # upper-diagonal
    dy = np.zeros((ny)) # right-hand side

    expl_term_temp_x = np.zeros((nx,ny))
    expl_term_temp_y = np.zeros((nx,ny))

    expl_term_old_x = np.zeros((nx,ny))
    expl_term_old_y = np.zeros((nx,ny))

    k_rel = relative_coeff(simtime,lambda_wb,k_wb)
    f_old = nonlinear_function(h_old,delta_x,delta_y,S_c,enne,k*k_rel,max_nlc)

    for inner_iter in range(0,max_inner_iter-1):

        if verbose_level == 2:

            print("Outer loop: dt = ",delta_t)

        elif verbose_level > 2:

            print("Outer loop: dt = ",delta_t)

        f = nonlinear_function(h_temp,delta_x,delta_y,S_c,enne,k*k_rel,max_nlc)

        # Evaluate the explicit term (the right-hand side of the linear system
        # for the diffusion in the x-direction)

        flux_temp_east,flux_temp_west,flux_temp_north,flux_temp_south = \
                            numerical_fluxes( h_temp , f , delta_x , delta_y )

        expl_term_temp_x[0,:] = flux_temp_east[0,:] / delta_x

        expl_term_temp_x[1:nx-1,:] = ( flux_temp_east[1:nx-1,:] \
                            - flux_temp_west[1:nx-1,:] ) / delta_x

        expl_term_temp_x[nx-1,:] = - flux_temp_west[nx-1,:] / delta_x

        expl_term_temp_y[0:nx,0] = flux_temp_north[0:nx,0] / delta_y

        expl_term_temp_y[0:nx,1:ny-1] = ( flux_temp_north[0:nx,1:ny-1] \
                            - flux_temp_south[0:nx,1:ny-1] ) / delta_y

        expl_term_temp_y[0:nx,ny-1] = - flux_temp_south[0:nx,ny-1] / delta_y

        #-----------------------------

        flux_old_east,flux_old_west,flux_old_north,flux_old_south = \
                            numerical_fluxes( h_old , f_old , delta_x , delta_y )

        expl_term_old_x[0,:] = flux_old_east[0,:] / delta_x

        expl_term_old_x[1:nx-1,:] = ( flux_old_east[1:nx-1,:] \
                            - flux_old_west[1:nx-1,:] ) / delta_x

        expl_term_old_x[nx-1,:] = - flux_old_west[nx-1,:] / delta_x

        expl_term_old_y[0:nx,0] = flux_old_north[0:nx,0] / delta_y

        expl_term_old_y[0:nx,1:ny-1] = ( flux_old_north[0:nx,1:ny-1] \
                            - flux_old_south[0:nx,1:ny-1] ) / delta_y

        expl_term_old_y[0:nx,ny-1] = - flux_old_south[0:nx,ny-1] / delta_y


        expl_term = - ( h_temp - h_old ) \
                    + delta_t * 0.5 * ( expl_term_temp_x + expl_term_temp_y ) \
                    + delta_t * 0.5 * ( expl_term_old_x + expl_term_old_y )   \
                    + A_c * delta_t


        # Solve the linearized system in the x-direction with the Thomas
        # algorithm.


        if ( Dirichlet[0] ) or ( Transient[0] ):

            bx[0] = 1.0;   # these coefficients are defined in order to impose the
            cx[0] = 0.0;   # Dirichlet boundary condition
            dx[0] = 0.0;   # boundary condition


        if ( Dirichlet[1] ) or ( Transient[1] ):

            first_row = 1

        else:

            first_row = 0

        if ( Dirichlet[2] ) or ( Transient[2] ):

            ax[nx-1] = 0.0   # these coefficients are defined in order to impose
            bx[nx-1] = 1.0   # the Dirichlet boundary condition
            dx[nx-1] = 0.0   # boundary condition

        if ( Dirichlet[3] ) or ( Transient[3] ):

            last_row = ny-1

        else:

            last_row = ny


        for j in range(first_row,last_row):

            if ( Neumann[0] ):

                # coefficients for the Neumann boundary condition at node (0,j)

                bx[0] = 1.0 + 0.25 * delta_t / delta_x**2 * \
                        ( f[0,j] + f[1,j] )

                cx[0] = - 0.25 * delta_t / delta_x**2 * ( f[1,j] + f[0,j] )
                dx[0] = expl_term[0,j]

            ax[1:nx-1] = - 0.25 * delta_t / delta_x**2 * ( f[0:nx-2,j] + f[1:nx-1,j] )

            bx[1:nx-1] = 1.0 + 0.25 * delta_t / delta_x**2 * \
               ( f[2:nx,j] + 2.0 * f[1:nx-1,j] + f[0:nx-2,j] )

            cx[1:nx-1] = - 0.25 * delta_t / delta_x**2 * ( f[2:nx,j] + f[1:nx-1,j] )

            dx[1:nx-1] = expl_term[1:nx-1,j]

            if ( Neumann[2] ):

               # coefficients for the Neumann boundary condition at node (nx,j)
               ax[nx-1] = - 0.25 * delta_t / delta_x**2 * ( f[nx-2,j] + f[nx-1,j] )

               bx[nx-1] = 1.0 + 0.25 * delta_t / delta_x**2 * \
                        ( f[nx-1,j] + f[nx-2,j] )

               dx[nx-1] = expl_term[nx-1,j]

            # call the tridiagonal solver
            # w[0:nx,j] = TDMAsolver(ax,bx,cx,dx)
            w[0:nx,j] = NEWsolver(ax,bx,cx,dx)


        # Solve in the y-direction

        if ( Dirichlet[0] ):

            first_column = 1
            h_temp[0,0:ny] = h_old[0,0:ny]

        elif ( Transient[0] ):

            first_column = 1
            h_temp[0,0:ny] = h_old[0,0:ny] + grow_rate[0]*delta_t

        elif ( Neumann[0] ):

            first_column = 0


        if ( Dirichlet[1] ):

            by[0] = 1.0
            cy[0] = 0.0
            dy[0] = 0.0

        elif ( Transient[1] ):

            by[0] = 1.0
            cy[0] = 0.0
            dy[0] = grow_rate(2)*delta_t

            h_temp[0:nx,0] = h_old[0:nx,0] + grow_rate[1]*delta_t


        if ( Dirichlet[2] ):

            last_column = nx-1
            h_temp[nx-1,0:ny] = h_old[nx-1,0:ny]

        elif ( Transient[2] ):

            last_column = nx-1
            h_temp[nx-1,0:ny] = h_old[nx-1,0:ny] + grow_rate[2]*delta_t

        elif ( Neumann[2] ):

            last_column = nx


        if ( Dirichlet[3] ):

            ay[ny-1] = 0.0
            by[ny-1] = 1.0
            dy[ny-1] = 0.0

        elif ( Transient[3] ):

            ay[ny-1] = 0.0
            by[ny-1] = 1.0
            dy[ny-1] = grow_rate[3]*delta_t

            h_temp[0:nx,ny-1] = h_old[0:nx,ny-1] + grow_rate[3]*delta_t;


        for i in range(first_column,last_column):

            if ( Neumann[1] ):

                # coefficients for the Neumann boundary condition at node [i,0]

                by[0] = 1.0 + 0.25 * delta_t / delta_y**2 * \
                        ( f[i,0] + f[i,1] )
                cy[0] = - 0.25 * delta_t / delta_y**2 * ( f[i,1] + f[i,0] )
                dy[0] = w[i,0]


            ay[1:ny-1] = - 0.25 * delta_t / delta_y**2 * ( f[i,0:ny-2] + f[i,1:ny-1] )
            by[1:ny-1] = 1.0 + 0.25 * delta_t / delta_y**2 * \
                         ( f[i,2:ny] + 2.0 * f[i,1:ny-1] + f[i,0:ny-2] )

            cy[1:ny-1] = - 0.25 * delta_t / delta_y**2 * ( f[i,2:ny] + f[i,1:ny-1] )

            dy[1:ny-1] = w[i,1:ny-1]


            if ( Neumann[3] ):

                # coefficients for the Neumann boundary condition at node [i,ny-1
                ay[ny-1] = - 0.25 * delta_t / delta_y**2 * ( f[i,ny-2] + f[i,ny-1] )
                by[ny-1] = 1.0 + 0.25 * delta_t / delta_y**2 * \
                           ( f[i,ny-1] + f[i,ny-2] )

                dy[ny-1] = w[i,ny-1]

            # call the tridiagonal solver
            # v[i,0:ny] = TDMAsolver(ay,by,cy,dy)
            v[i,0:ny] = NEWsolver(ay,by,cy,dy)

        residual = np.maximum(np.abs(v.min()), np.abs(v.max()) )

        if ( verbose_level >= 3 ):

            print('Inner loop. Iter = ' + str(i) + ' residual = ' + str(residual))


        h_temp[first_column:last_column,first_row:last_row] +=   \
                    v[first_column:last_column,first_row:last_row]

        if ( residual < res ):

            break


    h_new[0:nx,0:ny] = 1.0 * h_temp[0:nx,0:ny]


    return [h_new[0:nx,0:ny],residual]

def NEWsolver(a,b,c,d):

    import numpy as np
    import scipy.linalg as la

    # Create arrays and set values
    m = b.size
    ab = np.zeros((3,m))
    ab[0] = a
    ab[1] = b
    ab[2] = c

    return la.solve_banded ((1,1),ab,d)

def TDMAsolver(a,b,c,d):

    # TDMAsolver is tridiagonal matrix solver using the Thomas algorithm
    # (named after Llewellyn Thomas). It is a simplified form of Gaussian
    # elimination that can be used to solve tridiagonal systems of equations.
    # A tridiagonal system for n unknowns may be written as
    #
    #    a_i*x_{i-1} + b_i*x_i + c_i*x_{i+1} = d_i ,
    #
    # where a a_1=0 and c_n=0.
    # a, b, c, and d are the column 1xn arrays for the compressed tridiagonal
    # matrix.

    import numpy as np

    n = b.size  # n is the number of rows

    # Modify the first-row coefficients
    c[0] /= b[0]      # Division by zero risk.
    d[0] /= b[0]      # Division by zero would imply a singular matrix.

    x = np.zeros((n))

    for i in range(1,n):

        id = 1.0 / (b[i] - c[i-1] * a[i])  # Division by zero risk.
        c[i] *= id                         # Last value calculated is redundant.
        d[i] = (d[i] - d[i-1] * a[i]) * id

    # Now back substitute.
    x[n-1] = d[n-1]

    for i in reversed(range(n-1)):

        x[i] = d[i] - c[i] * x[i + 1]

    return x
