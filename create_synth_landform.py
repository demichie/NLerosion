def create_synth_landform(x_min, delta_x, x_max, y_min, delta_y, y_max,r1,r2,h1,h2,h3,a,b,symmetry):

    import numpy as np

    nx = np.int(np.ceil( ( x_max - x_min ) / delta_x ) ) + 1
    ny = np.int(np.ceil( ( y_max - y_min ) / delta_y ) ) + 1

    x = np.linspace(x_min,x_max,nx)
    y = np.linspace(y_min,y_max,ny)


    X,Y = np.meshgrid(x,y)


    h = np.zeros((ny,nx)) # initialization of the array of the altitudes

    for i in range(0,nx):

        for j in range(0,ny):

            if ( symmetry == 'radial' ):

                r = np.sqrt(x[i]**2/a**2+y[j]**2/b**2)

            elif ( symmetry == 'y'):

                r = x[i]/a

            elif ( symmetry == 'x'):

                r = y(j)/b

            else:

                sys.exit("Please define symmetry: radial, x, or y)")

            if ( r < r1 ):

                h[j,i] = h1

            elif ( r < r2 ):

                h[j,i] = h1 + ( r - r1 ) / (r2 - r1) * (h2 - h1)

            elif ( r < 1 ):

                h[j,i] = h2 + ( r - r2 ) / (1.0 - r2) * (h3 - h2)

            else:

                h[j,i] = h3


    # diffusion coefficient

    # there are 5 regions defined:
    # r<r1      k=k1
    # r1<r<r2   k changes linearly between k1 and k2
    # r2<r<r3   k=k2
    # r3<r<r4   k changes linearly between k2 and k1
    # r2>r4     k=k1

    r1 = 0.25 #relative radius
    r2 = 0.4
    r3 = 0.6
    r4 = 0.75

    k1 = 1
    k2 = 1 #one for rim

    k = np.zeros((ny,nx))  # initialization of the array of the erosion coefficients

    for i in range(0,nx):# a loop to define erosion coefficients for all the cells/columns of the grid

        for j in range(0,ny):

            r = np.sqrt(x[i]**2/a**2+y[j]**2/b**2)

            if ( r < r1 ):

                k[j,i] = k1

            elif ( r < r2 ):

                k[j,i] = k1 + ( r - r1 ) / (r2 - r1) * (k2 - k1) # j ( or first index) for the row, i for the columns

            elif ( r < r3 ):

                k[j,i] = k2;

            elif ( r < r4 ):

                k[j,i] = k2 + ( r - r3 ) / (r4 - r3) * (k1 - k2)

            else:

                k[j,i] = k1


    return (X, Y,h,k)
