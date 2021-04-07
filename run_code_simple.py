import numpy as np
from cinder_cone_and_k import cinder_cone_and_k
from ErosioNL import ErosioNL
from read_asc import read_asc

# Script to run the non linear diffusion code on various initial forms
# 2D non linear diffusion code is by 
# Mattia de' Michieli Vitturi (Istituto Nazionale di Geofisica e Vulcanologia)
# April 20, 2016

run_name = 'synthetic_cone'

# First, build an initial form

# Synthetic cinder cone (Hooper & Sheridan):

x_min = -500.0
delta_x = 10.0
x_max = 500.0
y_min = x_min
delta_y = delta_x
y_max = x_max 
[X, Y, h, k] = cinder_cone_and_k(x_min, delta_x, x_max, y_min, delta_y, y_max)

# Topography from ascii raster file
# ascii_file = 'synthetic_cone_init.asc'
# [X, Y, h, k,x_min,x_max,delta_x,y_min,y_max,delta_y] = read_asc(ascii_file)


# Save initial topography on ascii raster file
header = "ncols     %s\n" % h.shape[1]
header += "nrows    %s\n" % h.shape[0]
header += "xllcenter " + str(x_min) +"\n"
header += "yllcenter " + str(y_min) +"\n"
header += "cellsize " + str(delta_x) +"\n"
header += "NODATA_value -9999\n"

output_full = run_name + '_init.asc'

np.savetxt(output_full, np.flipud(h), header=header, fmt='%1.5f',comments='')
print(output_full+' saved')


# Initialize the mask defining the region to uplift/depression and/or tilt 

# no uplift and tilt

nx = h.shape[0]
ny = h.shape[1]


mask = np.zeros((nx,ny))


x0 = 0.0    # x-center of tilt
y0 = 0.0    # y-center of tilt
alfa = 0.0  # angle defining the tilt-axis
c0 = 0.0    # depression in the ragion defined by the mask 
c1 = 0.0    # no tilting along the axis parallel to the line alfa=0  
c2 = 0.0    # no tilting along the axis orthogonal to the line alfa=0 


# Now set up model input
# X, Y, h have been built already

final_time = 100.0   # final time in kilo years

delta_t_max = 10.00  # maximum time step
delta_t0 = 0.01      # initial time step

cr_angle = 33.0      # critical slope in degrees

enne = 2             # exponent for the nonlinearity of the model
                     # enne = inf  gives a linear model

k = 1.0              # m^2/kyr; the diffusion coefficient, a scalar or a (nx,ny) array
                     # the time unit is the same of the parameter 'final_time'

max_nlc = 10.0       # maximum value of the nonlinear coefficient

max_inner_iter = 100 # the maximum number of iteration for the inner loop. 
                     # <20 should be good;

res = 1.e-4          # [m] the residual required for the convergenge of the inner 
                     # loop.

bc = 'NNNN'  # the boundary condition : 'N' (Neumann) or 'D' (Dirichlet)
             # or 'T' (Transient). The order is S,W,N,E.
             # Dirichlet is fixed values (elevation in this case) 
             # Neumann is fixed gradient (flux null in this case)
             # Transient is elevation at the boundaries changing at
             # fixed rate, given by the parameter 'grow_rates'
    
gr = -100.0             
grow_rates = [ gr, gr, gr, gr]  # rate of change at the boundaries
                             # Used only when the b.c. is 'T'

n_output = 10  # number of output plotted



vx = X - x0
vy = Y - y0

alfarad = np.deg2rad(alfa)

A_c0 = c0 * mask
A_c1 = c1 * mask * ( vx * np.cos(alfarad) + vy * np.sin(alfarad) )
A_c2 = c2 * mask * ( vx * np.sin(alfarad) - vy * np.cos(alfarad) )
A_c = A_c0 + A_c1 + A_c2


verbose_level = 0   # level of output on screen (>1 for debug purposes)
plot_output_flag = 0
save_output_flag = 1

( h_new , h_diff , time_iter , max_angle , mean_angle  ) =       \
    ErosioNL(X,Y,h,final_time,delta_t_max,delta_t0,              \
    cr_angle,enne,k,max_nlc,max_inner_iter, res, bc, grow_rates, \
    run_name , n_output , A_c , verbose_level,save_output_flag,  \
     plot_output_flag)

