def nonlinear_function(h,delta_x,delta_y,S_c,enne,k,max_nlc):

    # NONLINEAR_FUNCTION This function evaluate the nonlinear function
    #
    # f(h) = k * ( 1 + ( max_nlc -1 ) * ( |grad(h)| / S_c )^n )
    #
    # k is the diffusivity constant, S_c is the tangent of the critical slope.
    # and max_nlc is the maximum value of the non-linear coefficient that
    # multiply k.

    import numpy as np

    nx = h.shape[0]
    ny = h.shape[1]

    h_x,h_y = np.gradient(h,delta_x,delta_y)

    grad_h = np.sqrt( h_x**2 + h_y**2 )

    ratio = np.minimum(np.ones((nx,ny)),grad_h/S_c) #gradient computed in each sells, so wnp.ones is an array

    if ( np.isinf(enne) ):

        non_linear_coeff = 1.0

    else:

        if ( enne == 0 ): #when n=0 we have linear flux(diffusion). enne = n!!!

            non_linear_coeff = max_nlc

        else:

            ratio2 = ( ( max_nlc - 1.0 ) / max_nlc ) ** ( 1.0 / enne ) * ratio #max_nlc=cmax
            non_linear_coeff = 1.0 / ( 1.0 - ratio2**enne )

    effe = k * non_linear_coeff

    return effe
