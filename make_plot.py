def make_plot(X,Y,Z,Z_diff,ax):

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import matplotlib.pyplot as plt
    import numpy as np

    ax.clear()

    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, 
                    cmap=cm.bone,linewidth=0, antialiased=False)

    if ( np.min(Z_diff) == np.max(Z_diff) ):

        Z_diff[0,0] = 0.001

    cset1 = ax.contour(X, Y, Z_diff, zdir='z', offset=-100, cmap=cm.coolwarm)
    cset2 = ax.contour(X, Y, Z, zdir='x', offset=-600, cmap=cm.coolwarm)
    cset3 = ax.contour(X, Y, Z, zdir='y', offset=600, cmap=cm.coolwarm)

    ax.set_xlabel('X')
    ax.set_xlim(-600, 600)
    ax.set_ylabel('Y')
    ax.set_ylim(-600, 600)
    ax.set_zlabel('Z')
    ax.set_zlim(-50, 200)


    plt.show()
    plt.pause(0.01)
