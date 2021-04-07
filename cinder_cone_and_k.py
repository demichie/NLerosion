def cinder_cone_and_k(x_min, delta_x, x_max, y_min, delta_y, y_max):

    import numpy as np

    nx = np.int(np.ceil( ( x_max - x_min ) / delta_x ) ) + 1 
    ny = np.int(np.ceil( ( y_max - y_min ) / delta_y ) ) + 1

    x = np.linspace(x_min,x_max,nx)
    y = np.linspace(y_min,y_max,ny)


    X,Y = np.meshgrid(x,y)

    # parameters defining the cone

    # r is the distance from the center of the cone
    # the top is at r2 
    # there are 4 regions defined:
    # r<r1      h = h1
    # r1<r<r2   h varying linearly between h1 and h2
    # r2<r<r3   h varying linearly between h2 and h3
    # r<r3      h = h3

    r1 = 50
    r2 = 100
    r3 = 350

    h1 = 130
    h2 = 160
    h3 = 0

    h = np.zeros((ny,nx)) # initialization of the array of the altitudes

    for i in range(0,nx):
    
        for j in range(0,ny):
        
            r = np.sqrt( x[i]**2 + y[j]**2 )
        
            if ( r < r1 ):
            
                h[j,i] = h1
            
            elif ( r < r2 ):
            
                h[j,i] = h1 + ( r - r1 ) / (r2 - r1) * (h2 - h1)
            
            elif ( r < 350.0 ):
            
                h[j,i] = h2 + ( r - r2 ) / (r3 - r2) * (h3 - h2)
            
            else:
            
                h[j,i] = h3
            

    # diffusion coefficient

    # there are 5 regions defined:
    # r<r1      k=k1
    # r1<r<r2   k changes linearly between k1 and k2
    # r2<r<r3   k=k2
    # r3<r<r4   k changes linearly between k2 and k1
    # r2>r4     k=k1
    
    r1 = 40
    r2 = 50
    r3 = 150
    r4 = 160
    
    k1 = 10
    k2 = 10
    
    k = np.zeros((ny,nx))  # initialization of the array of the erosion coefficients

    for i in range(0,nx):
    
        for j in range(0,ny):
        
            r = np.sqrt( x[i]**2 + y[j]**2 )
        
            if ( r < r1 ):
            
                k[j,i] = k1
            
            elif ( r < r2 ):
            
                k[j,i] = k1 + ( r - r1 ) / (r2 - r1) * (k2 - k1)
            
            elif ( r < r3 ):
            
                k[j,i] = k2;
                        
            elif ( r < r4 ):
            
                k[j,i] = k2 + ( r - r3 ) / (r4 - r3) * (k1 - k2)
            
            else:
            
                k[j,i] = k1

    
    return (X, Y,h,k)


