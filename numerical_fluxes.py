def numerical_fluxes( h , f , delta_x , delta_y ):

    # NUMERICAL_FLUXES This function avaluates the fluxes at the boundaries
    # of each cell (east,west,south,north) with a finite-difference scheme.

    import numpy as np
 
    nx = h.shape[0]
    ny = h.shape[1]


    flux_east = np.zeros((nx,ny))
    flux_west = np.zeros((nx,ny))
    flux_north = np.zeros((nx,ny))
    flux_south = np.zeros((nx,ny))


    flux_west[1:nx,:] = 0.5 * ( f[0:nx-1,:] + f[1:nx,:] ) * \
    ( h[1:nx,:] - h[0:nx-1,:] ) / delta_x

    flux_east[0:nx-1,:] = flux_west[1:nx,:]

    flux_south[:,1:ny] = 0.5 * ( f[:,0:ny-1] + f[:,1:ny] ) * \
    ( h[:,1:ny] - h[:,0:ny-1] ) / delta_y

    flux_north[:,0:ny-1] = flux_south[:,1:ny]

    return (flux_east,flux_west,flux_north,flux_south)

