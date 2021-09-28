def read_asc(ascii_file):

    import numpy as np
    from linecache import getline

    source1 = ascii_file
    # Parse the header using a loop and
    # the built-in linecache module
    hdr = [getline(source1, i) for i in range(1,7)]
    values = [float(h.strip('\n').strip().split(" ")[-1]) \
     for h in hdr]
    cols,rows,xll,yll,cell,nd = values
    cols = int(cols)
    rows = int(rows)
    delta_x = cell
    delta_y = cell

    values = [(h.split(" ")[0]) for h in hdr]
    s1,s2,s3,s4,s5,s6 = values

    if ( s3=='xllcorner'):

        x_min = xll+0.5*delta_x
        y_min = yll+0.5*delta_y
    
    elif ( s3=='xllcenter'):
    
        x_min = xll
        y_min = yll

    x_max = x_min+(cols-1)*delta_x
    x = np.linspace(x_min,x_max,cols)

    y_max = y_min+(rows-1)*delta_y
    y = np.linspace(y_min,y_max,rows)

    X,Y = np.meshgrid(x,y)

    # Load the dem into a numpy array
    h = np.flipud(np.loadtxt(source1, skiprows=6))

    h[h==nd]=np.nan

    return (X, Y, h, x_min, x_max, delta_x, y_min, y_max, delta_y)


