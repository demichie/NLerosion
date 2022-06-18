import numpy as np
    
def regrid(xin, yin, fin, xl, xr , yl, yr):

    """
    Interpolation from a regular grid to a second regular grid 

    @params:
        xin      - Required : original grid X values (1D Dble) 
        yin      - Required : original grid Y values (1D Dble) 
        fin      - Required : original grid Z values (2D Dble) 
        xout     - Required : new grid X values (2D Dble) 
        yout     - Required : new grid Y values (2D Dble) 
    """   
         
    nXin = len(xin)-1
    nYin = len(yin)-1

    dXin = xin[1] - xin[0]
    dYin = yin[1] - yin[0]
    
    ix1 = np.maximum( 0 , np.ceil(( xl - xin[0] ) / dXin ).astype(int) -1 )
    ix2 = np.minimum( nXin , np.ceil( ( xr -xin[0] ) / dXin ).astype(int) )
        
    iy1 = np.maximum( 0 , np.ceil( ( yl - yin[0] ) / dYin ).astype(int) -1 )
    iy2 = np.minimum( nYin , np.ceil( ( yr - yin[0] ) / dYin ).astype(int) )

    fout = 0.0

    for ix in range(ix1,ix2):
        
       alfa_x = ( np.minimum(xr,xin[ix+1]) - np.maximum(xl,xin[ix]) ) / ( xr - xl )

       for iy in range(iy1,iy2):
                
          alfa_y = ( np.minimum(yr,yin[iy+1]) - np.maximum(yl,yin[iy]) ) / ( yr - yl )
          
          #print('1',fout == -9999)
          #print('alfa_x',alfa_x)
          #print('yl,yr',yl,yr)
          #print('iy',iy)
          #print('yin[iy+1]',yin[iy+1])
          #print('alfa_y',alfa_y)
          #print('2',alfa_x * alfa_y > 0.0)
          #print('3',fin[iy,ix] == - 9999)
          if (fout == -9999) or ( ( alfa_x * alfa_y > 0.0 ) and ( fin[iy,ix] == - 9999 ) ):
          
              fout = -9999
              
          else:    
          
              # print('alfa_x*alfa_y',alfa_x*alfa_y)
              fout = fout + alfa_x * alfa_y * fin[iy,ix]

    return fout
 

def regrid2Dgrids(xin,yin,Zin,Xout,Yout):
    """
    Interpolation from a regular grid to a second regular grid 

    @params:
        xin      - Required : original grid-centers X values (2D Dble) 
        yin      - Required : original grid-centers Y values (2D Dble) 
        Zin      - Required : original grid-centers Z values (2D Dble) 
        xout     - Required : new grid-centers X values (2D Dble) 
        yout     - Required : new grid-centers Y values (2D Dble) 
    """

    cellin = xin[0,1]-xin[0,0]
    
    # print('cellin',cellin)

    x1 = [x-0.5*cellin for x in xin[0,:]]
    x1.append(x1[-1]+cellin)
        
    y1 = [y-0.5*cellin for y in yin[:,0]]
    y1.append(y1[-1]+cellin)

    #x1 = np.array(x1)
    #y1 = np.array(y1)

    # print(x1)
    # print(y1)


    nrows = Xout.shape[0]
    ncols = Xout.shape[1]
        
    cellout = Xout[0,1] - Xout[0,0]
    lxout = Xout[0,0] - 0.5*cellout
    lyout = Yout[0,0] - 0.5*cellout

    Zout = np.zeros_like(Xout)

    for j in range(0,ncols):
           
        xl = lxout + (j)*cellout
        xr = lxout + (j+1)*cellout
              
        for k in range(0,nrows):
             
            yl = lyout + (k)*cellout
            yr = lyout + (k+1)*cellout
             
            Zout[k,j] = regrid( x1 , y1 , Zin , xl , xr , yl , yr )
        
    return Zout    

def read_asc(ascii_file,cellsize='none'):

    import numpy as np
    from linecache import getline

    print('Reading file: ',ascii_file)
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

    print('Reading file completed')
    
    if cellsize == 'none':
    
        cellsize = cell
    
    if cellsize != cell:
    
        print('Resampling DEM')
    
        xl = x_min - 0.5*cell
        xr = x_max + 0.5*cell

        yl = y_min - 0.5*cell
        yr = y_max + 0.5*cell

        # edges of resampled grid   
        xnew = np.arange(xl,xr,cellsize)
        ynew = np.arange(yl,yr,cellsize)
        
        # centers of resamples grid
        Xnew,Ynew = np.meshgrid(xnew[:-1]+0.5*cellsize,ynew[:-1]+0.5*cellsize)
        
        hnew = regrid2Dgrids(X,Y,h,Xnew,Ynew)
        print('Resampling completed')
            
        X = Xnew
        Y = Ynew
        h = hnew
        x_min = np.amin(X)
        x_max = np.amax(X)
        y_min = np.amin(Y)
        y_max = np.amax(Y)
        
        delta_x = cellsize
        delta_y = cellsize
        

    return (X, Y, h, x_min, x_max, delta_x, y_min, y_max, delta_y)


