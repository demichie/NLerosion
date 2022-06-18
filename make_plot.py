def make_plot(X,Y,Z,Z_init,h_min,h_max,simtime,cr_angle,run_name,iter,
              plot_show_flag,flank_mask,dist):

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib import colors
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import LightSource
    import numpy.ma as ma
    import pandas as pd
    
    import imp

    try:
        imp.find_module('mayavi')
        found = True
    except ImportError:
        found = False
        
    if found:
    
        from mayavi import mlab    

    plt.ion()
    plt.rcParams.update({'font.size':8})

     
    plt.close('all') 
     
    Z_diff = Z-Z_init 
        
    fig = plt.figure()
    fig.set_size_inches(11,7)
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    plt.tight_layout(pad=4, w_pad=4, h_pad=4)
    time_text = ax3.text(0,0, 'time ='+"{:8.2f}".format(simtime)+'s')
    
    delta_x = X[0,1]-X[0,0]
    delta_y = Y[1,0]-Y[0,0]
    
    x_min = X[0,0] - 0.5*delta_x
    x_max = X[0,-1] + 0.5*delta_x
    
    y_min = Y[0,0] - 0.5*delta_y
    y_max = Y[-1,0] + 0.5*delta_y

    Z_x,Z_y = np.gradient(Z,delta_x,delta_y)

    grad_Z = np.sqrt( Z_x**2 + Z_y**2 )

    slope = ( np.arctan(grad_Z)*180.0/np.pi )
    
    my_col = cm.jet(slope/cr_angle)

    norm = colors.Normalize(vmin=0.0, vmax=cr_angle)

    
    ls = LightSource(azdeg=315, altdeg=45)

    extent = [x_min, x_max, y_min, y_max]
    ax1.imshow(ls.hillshade(Z,vert_exag=1.0,
                 dx=delta_x,dy=delta_y),cmap='gray',extent=extent,
                 origin='lower') 

    if(np.isnan(flank_mask).any()):
    
        flank_plot = np.zeros_like(flank_mask)
        flank_plot[flank_mask>0] = 1.0

        ax1.imshow(flank_plot,cmap='gray',extent=extent,alpha=0.2,
                 origin='lower') 
                                        
    nx = X.shape[1]
    ny = X.shape[0]  
    
    idx1 = int(np.floor(nx/4)) 
    idx2 = int(np.floor(nx/2)) 
    idx3 = int(np.floor(3*nx/4)) 
    
    dh = h_max-h_min
    
    ax1.plot(X[idx1,:],Y[idx1,:], 'blue')
    ax1.plot(X[idx2,:],Y[idx2,:], 'red')
    ax1.plot(X[idx3,:],Y[idx3,:], 'green')

    if ( np.min(Z_diff) == np.max(Z_diff) ):
    
        Z_diff[0,0] = 0.001
        Z_diff[-1,-1] = -0.001
        
   
    ax1.set_xlabel('x [m]')
    ax1.set_xlim(x_min,x_max)
    ax1.set_ylabel('y [m]')
    ax1.set_ylim(y_min,y_max)
    
    z_range = np.amax(np.abs(Z_diff))   
    cnt = ax2.imshow(Z_diff, cmap='seismic', extent=extent, vmin=-z_range, 
                     vmax=z_range, origin='lower')
    
    clb = plt.colorbar(cnt,ax=ax2)
    clb.set_label('Delta h [m]')
    # colorbar(cnt);
    
    ax2.set_xlabel('x [m]')
    ax2.set_xlim(np.amin(X),np.amax(X))
    ax2.set_ylabel('y [m]')
    ax2.set_ylim(np.amin(Y),np.amax(Y))
    

    
    l1, = ax3.plot(X[idx1,:],Z[idx1,:], 'b-')
    l2, = ax3.plot(X[idx2,:],Z[idx2,:], 'r-')
    l3, = ax3.plot(X[idx3,:],Z[idx3,:], 'g-')
    
    ax3.legend((l1, l2, l3), ("y="+str(Y[idx1,0])+' m',"y="+str(Y[idx2,0])+' m', \
                             "y="+str(Y[idx3,0])+' m'), loc='upper right', shadow=True)

    
    ax3.plot(X[idx1,:],Z_init[idx1,:],'b--')
    ax3.plot(X[idx2,:],Z_init[idx2,:],'r--')
    ax3.plot(X[idx3,:],Z_init[idx3,:],'g--')
           
    ax3.set_xlabel('x [m]')      
    ax3.set_ylabel('z [m]')    
    
    x_min, x_max = ax3.get_xlim()
        
    y_min, y_max = ax3.get_ylim()
    
    time_text.set_position((x_min+0.05*(x_max-x_min), y_min+0.9*(y_max-y_min)))
    time_text.set_text('time ='+"{:8.2f}".format(simtime)+'s')

    flank_dist = np.zeros_like(flank_mask)
    flank_dist[flank_mask>0] = dist[flank_mask>0]

    masked_slope = ma.masked_array(slope,flank_dist == 0)
    slope_flat = masked_slope[masked_slope.mask == False].flatten()
     
    masked_dist = ma.masked_array(flank_dist, flank_dist == 0)
    dist_flat = masked_dist[masked_dist.mask == False].flatten()

    df = pd.DataFrame({'dist': dist_flat, 'slope': slope_flat})

    n_dist = 10
    # create the bins (intervals) for the sectors based on the distance
    dist_bins = np.linspace(np.min(dist_flat), np.max(dist_flat),
                                n_dist + 1)
    dist_half = 0.5 * (dist_bins[0:-1] + dist_bins[1:])

    # print('dist_bins', dist_bins)

    # group the dataframe elements in sub-domains by using partition
    groups = df.groupby([pd.cut(df.dist, dist_bins)])

    # compute the mean slope in each sub-domain. It is a 2D numpy array
    slp_mean = np.array(groups['slope'].mean())
    slp_std = np.array(groups['slope'].std())

    # print('slp_mean', slp_mean)
    # print('slp_std', slp_std)
        
    ax4.errorbar(dist_half, slp_mean, yerr=slp_std)
    ax4.set_xlabel('distance from top [m]')
    ax4.set_ylabel('slope [degrees]')
    ax4.set_ylim(0,35)
                
 
    """
        
    # Get the histogramp
    cm = plt.cm.get_cmap('jet')

    Yh,Xh = np.histogram(slope_flat,density=True)
    x_span = Xh.max()-Xh.min()
    C = [cm(x/cr_angle) for x in Xh]

    ax4.bar(Xh[:-1]+0.5*(Xh[1]-Xh[0]),Yh,color=C,width=Xh[1]-Xh[0])      

    ax4b = ax4.twinx()
    
    # plot the cumulative histogram
    n_bins = 50
    ax4b.hist(slope_flat, n_bins, density=True, histtype='step',
                           cumulative=True, label='Empirical')

    ax4.set_xlabel('slope [degrees]')
    ax4.set_ylabel('probability density function')
    ax4b.set_ylabel('cumulative distribution function')

    """

    frame_name = run_name + '_{0:03}'.format(iter) + '.png'
    plt.savefig(frame_name,dpi=200)

    frame_name = run_name + '_{0:03}'.format(iter) + '.pdf'
    plt.savefig(frame_name)

    if plot_show_flag:

        plt.show()
        plt.pause(0.01)
    
    return time_text
