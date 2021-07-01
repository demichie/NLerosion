def ErosioNL( X, Y, h, final_time, delta_t_max,delta_t0, cr_angle, \
    enne, k, max_nlc, max_inner_iter, res , bc, grow_rate ,        \
    run_name , n_output , A_c , verbose_level, save_output_flag,   \
    plot_output_flag):


    # EROSIONL This function solve for the nonlinear diffusion problem 
    #
    # dh/dt = div( f(h) grad(h) )
    #
    # where
    #
    # f(h) = k * ( 1 + ( max_nlc -1 ) * ( |grad(h)| / S_c )^n )
    #
    # k is the diffusivity constant, S_c is the tangent of the critical slope.
    # and max_nlc is the maximum value of the non-linear coefficient that
    # multiply k.
    #
    # - X and Y are the (nx,ny) array of the coordinates as generated by the
    #   matlab meshgrid command;
    # - h is the (nx,ny) initial elevation array;
    # - final_time is the duration of the run;
    # - delta_t0 is the maximum time step;
    # - angle is the critical slope;
    # - enne is the exponent of the model;
    # - k is the diffusion goefficient, a scalar or a (nx,ny) array;
    # - max_nlc is the maximum value of the nonlienar coefficient;
    # - max_inner_iter is the maximum number of iteration for the inner loop;
    # - res is the residual required for the convergenge of the inner loop;
    # - bc is a 1x4 array with the boundary conditions at the S,W,N and E
    #   boundaries and the value can be: 'N' (Neumann), 'D' (Dirichlet) or
    #   'T' transient.
    # - grow_rate is a 1x4 array with the rates of change at the boundaries
    #   when the 'transient' boundary condition is applied;
    # - run_name is the name given to the output files
    # - n_output is the number of output plotted on screen.
    #
    # The model is a modification of the model presented in:
    #   Pelletier, J. # Cline, M., Nonlinear slope-dependent sediment transport
    #   in cinder cone evolution, Geology, Geological Society of America, 2007,
    #   35, 1067.
    # The numerical scheme adopted is described in:
    #   Tian, Z. & Ge, Y., A fourth-order compact ADI method for solving
    #   two-dimensional unsteady convection-diffusion problems, Journal of
    #   computational and applied mathematics, Elsevier, 2007, 198, 268-286.
    #
    # See also GRADIENT

    import numpy as np

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import matplotlib.pyplot as plt

    from advance_time import advance_time
    from nonlinear_function import nonlinear_function
    from make_plot import make_plot
    import netCDF4 
    import time
    import sys

    nx = h.shape[0]
    ny = h.shape[1]

    h_old = np.zeros((nx,ny))
    h_new = np.zeros((nx,ny))
    h_temp = np.zeros((nx,ny))
    h_temp_half = np.zeros((nx,ny))
    h_temp_full = np.zeros((nx,ny))

    if ( np.isscalar(k) ):
    
        k = k * np.ones((nx,ny))
    
    simtime = 0.0
    iter = 0

    iteration_draw_times = np.linspace(final_time/n_output,final_time,n_output)

    h_init = h

    if ( plot_output_flag ):

        plt.ion()
    
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        make_plot(X,Y,h,h-h_init,ax)

        frame_name = run_name + '_{0:03}'.format(0) + '.png'
        plt.savefig(frame_name,dpi=200)


    if ( save_output_flag):

        # Save initial topography on ascii raster file
        header = "ncols     %s\n" % h.shape[1]
        header += "nrows    %s\n" % h.shape[0]
        header += "xllcenter " + str(np.amin(X)) +"\n"
        header += "yllcenter " + str(np.amin(Y)) +"\n"
        header += "cellsize " + str(np.abs(X[1,2]-X[1,1])) +"\n"
        header += "NODATA_value -9999\n"

        output_full = run_name + '_{0:03}'.format(0) + '.asc'

        np.savetxt(output_full, np.flipud(h), header=header, fmt='%1.5f',comments='')


    # create netcdf4 file
    ncfilename = run_name+'.nc'

    ncfile = netCDF4.Dataset(ncfilename,mode='w',format='NETCDF4') 

    x_dim = ncfile.createDimension('x', nx) 
    y_dim = ncfile.createDimension('y', ny) 
    time_dim = ncfile.createDimension('time', None) # unlimited axis (can be appended to).

    ncfile.title = run_name+' output'

    ncfile.Conventions = "CF-1.0"
    ncfile.subtitle="My model data subtitle"
    ncfile.anything="write anything"


    x = ncfile.createVariable('x', np.float64, ('x',),zlib=True)
    x.long_name = 'x dim'
    x.units = 'meters'

    y = ncfile.createVariable('y', np.float64, ('y',),zlib=True)
    y.long_name = 'y dim'
    y.units = 'meters'

    t = ncfile.createVariable('time', np.float64, ('time',))
    t.long_name = 'Time'
    t.units = 'seconds'

    z = ncfile.createVariable('z',np.float64,('time','y','x'),zlib=True) # note: unlimited dimension is leftmost
    z.standard_name = 'flow thickness' # this is a CF standard name
    z.units = 'meters' 
    
    dz_dt = ncfile.createVariable('dz_dt',np.float64,('time','y','x'),zlib=True) # note: unlimited dimension is leftmost
    dz_dt.standard_name = 'erosion/deposit rate' # this is a CF standard name
    dz_dt.units = 'meters/seconds' 

    dz_tot = ncfile.createVariable('dz_tot',np.float64,('time','y','x'),zlib=True) # note: unlimited dimension is leftmost
    dz_tot.standard_name = 'flow erosion and deposit thickness' # this is a CF standard name
    dz_tot.units = 'meters' 

    t[iter] = simtime
    z[iter,:,:] = h
    dz_dt[iter,:,:] = h - h_init
    dz_tot[iter,:,:] = h - h_init

    x[:] = X[0,:]
    y[:] = Y[:,0]

    delta_x = np.abs(X[1,2]-X[1,1])
    delta_y = np.abs(Y[2,1]-Y[1,1])

    h_x,h_y = np.gradient(h,delta_x,delta_y)
    grad_h = np.sqrt( h_x**2 + h_y**2 )

    max_angle = np.zeros((n_output)) 
    mean_angle = np.zeros((n_output))
    time_iter = np.zeros((n_output))


    max_angle[iter] = np.arctan(np.max(np.max(grad_h)))*180/np.pi
    mean_angle[iter] = np.arctan(np.mean(grad_h))*180/np.pi

    time_iter[iter] = simtime

    print('Maximum/mean slope of the DEM = ' + str(max_angle[iter]))

    print('Mean slope of the DEM = ' + str(mean_angle[iter]))

    if ( cr_angle < max_angle[iter] ):
    
        print('Critical slope = ' + str(cr_angle))
        print('WARNING: critical slope < maximum slope of the DEM')
    

    S_c = np.tan(cr_angle/180.0*np.pi)

    # set the variables for the boundary conditions

    Neumann = np.zeros((4), dtype=bool)
    Dirichlet = np.zeros((4), dtype=bool)
    Transient = np.zeros((4), dtype=bool)


    for i in range(0,4):
    
        if ( bc[i] == 'N' ): 
        
            Neumann[i] = True
            Dirichlet[i] = False
            Transient[i] = False
            
        elif ( bc[i] == 'D' ):
            
            Dirichlet[i] = True
            Neumann[i] = False
            Transient[i] = False
            
        elif ( bc[i] == 'T' ):
            
            Transient[i] = True
            Dirichlet[i] = False
            Neumann[i] = False

        else:
            
            print('WARNING: wrong boundary condition, set to Neumann')
            Dirichlet[0:3] = False
            Neumann[0:3] = True
            



    f = nonlinear_function(h,delta_x,delta_y,S_c,enne,k,max_nlc)

    delta_t_stab = 0.5 * delta_x * delta_y / np.max(f)

    delta_t = np.minimum( delta_t0 , delta_t_stab )

    print('Delta t for stability = ', str(delta_t_stab))


    err = 1.0
    err_old = 1.0

    if sys.version_info >= (3, 0):

        start = time.process_time()

    else:

        start = time.clock()


    while ( simtime < final_time):
    
        delta_t_orig = delta_t
    
        coeff_t = np.maximum( 0.85 , np.minimum( 1.05,err**( -0.35 ) * err_old**( 0.2 ) ) )
    
        delta_t_orig = delta_t_orig * coeff_t
        delta_t_orig = min(delta_t_orig,iteration_draw_times[iter]-simtime)

        err = 1.0

        h_old[0:nx,0:ny] = h[0:nx,0:ny]
    
        while ( err >= 1.0 ):
                
            # integration with two steps with half step-size
        
            delta_t = 0.5 * delta_t_orig
        
            # first half step

        
            (h_temp[0:nx,0:ny],residual) = \
                    advance_time( h_old[0:nx,0:ny] , h_old[0:nx,0:ny] , delta_t , delta_x , \
                    delta_y , S_c , enne , k , max_nlc , max_inner_iter , res , \
                    Dirichlet, Transient , Neumann , grow_rate , A_c , verbose_level )

            # print '1st half ',np.min(h_temp),np.max(h_temp)
        
            # second half step
            
            h_guess = h_temp + ( h_temp - h_old)
        
            (h_new[0:nx,0:ny],residual) = \
                    advance_time( h_guess[0:nx,0:ny] , h_temp[0:nx,0:ny] , delta_t , delta_x , \
                    delta_y , S_c , enne , k , max_nlc , max_inner_iter , res , \
                    Dirichlet, Transient , Neumann , grow_rate , A_c , verbose_level )

            h_temp_half[0:nx,0:ny] = h_new[0:nx,0:ny]
        
            # print '1st half ',np.min(h_temp),np.max(h_temp)
            # print '2nd half ',np.min(h_temp_half),np.max(h_temp_half)

            # integration with one full step
        
            delta_t = delta_t_orig

            h_guess[0:nx,0:ny] = h_temp_half[0:nx,0:ny]

            (h_temp_full[0:nx,0:ny],residual) = \
                    advance_time( h_guess[0:nx,0:ny] , h_old[0:nx,0:ny] , delta_t , delta_x , \
                    delta_y , S_c , enne , k , max_nlc , max_inner_iter , res , \
                    Dirichlet, Transient , Neumann , grow_rate , A_c , verbose_level )
                
            # print '2nd half ',np.min(h_temp_half),np.max(h_temp_half)
            # print 'full step ',np.min(h_temp_full),np.max(h_temp_full)

            # check the error
        
            delta_t = delta_t_orig
            
            delta_h = np.abs( h_temp_half - h_temp_full )

            atol = res
            rtol = res

            err = np.sqrt( nx*ny * np.sum( ( delta_h / ( atol + np.abs(h) * rtol ) )**2 ) )
        
            if ( verbose_level >= 1 ):

                print('Accuracy = ' + str(err))
                
            if ( err >= 1.0):
            
                delta_t_orig = delta_t_orig * 0.975
                if ( verbose_level >= 1 ):
                    print('delta_t_orig =' ,delta_t_orig)

        err_old = err
        
        h_init  = h_init + A_c * delta_t
    
        h_diff = h_temp - h_init
    
        h[0:nx,0:ny] = h_temp_half[0:nx,0:ny]
    
        simtime = simtime + delta_t
    
        print('Time = %.5f delta_t = %.5f err = %.3f' % (simtime,delta_t,err))        

        if ( verbose_level >= 1):

            # gradient is computed using second order accurate central differences in the interior 
            # points and either first or second order accurate one-sides (forward or backwards) 
            # differences at the boundaries
            h_x,h_y = np.gradient(h,delta_x,delta_y)

            # magnitude of gradient
            grad_h = np.sqrt( h_x**2 + h_y**2 )

            max_angle[iter] = np.arctan(np.max(np.max(grad_h)))*180/np.pi
            mean_angle[iter] = np.arctan(np.mean(grad_h))*180/np.pi
            
            print('Maximum/mean slope of the DEM = ' + str(max_angle[iter]) \
                  + str(mean_angle[iter]))

        if ( simtime >= iteration_draw_times[iter] ):
            
            iter += 1
            print('Saving output '+str(iter))

            # add output to netcdf file
            t[iter] = iteration_draw_times[iter-1]
            z[iter,:,:] = h
            dz_dt[iter,:,:] = ( h - h_old ) / delta_t
            dz_tot[iter,:,:] = h - h_init
            
            if ( plot_output_flag):

                make_plot(X,Y,h_temp_half,h_temp_half-h_init,ax)
                            
                frame_name = run_name + '_{0:03}'.format(iter) + '.png'
                plt.savefig(frame_name,dpi=200)

            if ( save_output_flag):

                # Save topography on ascii raster file
                header = "ncols     %s\n" % h.shape[1]
                header += "nrows    %s\n" % h.shape[0]
                header += "xllcenter " + str(np.amin(X)) +"\n"
                header += "yllcenter " + str(np.amin(Y)) +"\n"
                header += "cellsize " + str(np.abs(X[1,2]-X[1,1])) +"\n"
                header += "NODATA_value -9999\n"

                output_full = run_name + '_{0:03}'.format(iter) + '.asc'
                print('Saving ascii rater file '+output_full)

                np.savetxt(output_full, np.flipud(h), header=header, fmt='%1.5f',comments='')


    if sys.version_info >= (3, 0):
    
        elapsed = (time.process_time() - start)

    else:
    
        elapsed = (time.clock() - start)

    print ('Time elapsed ' + str(elapsed) + ' sec.')
    print ('')
                 
    ncfile.close(); print('Dataset is closed!')


    h_new[0:nx,0:ny] = h[0:nx,0:ny]


    total_area = ( X[0,ny-1] - X[0,0] ) * ( Y[nx-1,0] - Y[0,0] )
    print('Total volume eroded (m^3) = ' + str(np.sum(h_new-h_init)*total_area))

    print('Max. deposition (m) = ' + str(np.maximum(0.0,np.max(h_diff))))

    print('Max. erosion (m) = ' + str(np.maximum(0.0,-np.min(h_diff))))

    plt.ioff()
    plt.show()


    return ( h_new , h_diff , time_iter , max_angle , mean_angle )







