#!/usr/bin/env python3

import shutil
import os
import numpy as np
from create_synth_landform import create_synth_landform
from ErosioNL import ErosioNL
from read_asc import read_asc

# Script to run the non linear diffusion code on various initial forms
# 2D non linear diffusion code is by
# Mattia de' Michieli Vitturi (Istituto Nazionale di Geofisica e Vulcanologia)
# April 20, 2016


# First, build an initial form


run_name = 'synth_cone'

# Synthetic cinder cone (Hooper & Sheridan):

# computational grid parameters
x_min = -500.0
delta_x = 10.0
x_max = 500.0
y_min = x_min
delta_y = delta_x
y_max = x_max

# parameters defining the cone

# semi-axis of the base of the cinder cone:
# a - semiaxis in the x-direction
# b - semiaxis in the y-direction
a = 400
b = 400

# r is defined as sqrt( (x/a)^2 + (y/b)^2 )
# r = 1 defines the base of the cone
# there are 4 regions defined:
# r<r1      h = h1                                (flat region inside the crater)
# r1<r<r2   h varying linearly between h1 and h2  (crater)
# r2<r<r3   h varying linearly between h2 and h3  (cone slope)
# r<1       h = h3                                (flat ragion outside the cone)

r1 = 0.25 #radius for elevation profile
r2 = 0.5

h1 = 100
h2 = 130
h3 = 0

symmetry = 'radial'

[X, Y, h, k] = create_synth_landform(x_min, delta_x, x_max, y_min, \
                                     delta_y, y_max,r1,r2,h1,h2,h3,a,b,symmetry)

"""
# Linear case:
run_name = 'synth_linear'

# computational grid parameters
x_min = -500.0
delta_x = 10.0
x_max = 500.0
y_min = x_min
delta_y = delta_x
y_max = x_max

# parameters defining the linear profile

# semi-axis of the base of the cinder cone:
# a - position of last slope change
# b - not used
a = 200
b = 200

# r is defined as x/a
# r = 1 defines the end of slope change zone
# there are 4 regions defined:
# r<=r1      h = h1                                (first flat region)
# r1<r<r2    h varying linearly between h1 and h2  (first constant slope region)
# r2<r<r1    h varying linearly between h2 and h3  (second constant slope region)
# r>=1       h = h3                                (final flat ragion)

r1 = 0.25
r2 = 0.5

h1 = 130
h2 = 160
h3 = 0

symmetry = 'y'

[X, Y, h, k] = synth_landform(x_min, delta_x, x_max, y_min, \
                              delta_y, y_max,r1,r2,h1,h2,h3,a,b,symmetry)
"""

# Topography from ascii raster file

"""
run_name = 'SPcrater'

ascii_file = 'sp_crater.asc'
[X, Y, h,x_min,x_max,delta_x,y_min,y_max,delta_y] = read_asc(ascii_file)
"""


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

final_time = 1000.0   # final time in kilo years

delta_t_max = 10.00  # maximum time step
delta_t0 = 0.01      # initial time step

cr_angle = 33.0      # critical slope in degrees

enne = 2.0           # exponent for the nonlinearity of the model
                     # enne = np.inf  gives a linear model

#k = 1.0              # m^2/kyr; the diffusion coefficient, a scalar or a (nx,ny) array
                     # the time unit is the same of the parameter 'final_time'

# The diffusion coefficient is multiplied by a relative coefficient k_rel, which
# is defined by a Weibull cumulative distribution with parameters lambda_wb and k_wb
   
lambda_wb = 20       # lambda_wb is the time scale parameter of the Weibull cumulative 
                     #   distribution (for t=lambda_wb, k_rel=(e-1)/e). If lambda_wb=0,
                     #   then k_rel is set to 1.0
                                          
k_wb = 5             # k_wb is the shape parameter of the Weibull distribution;

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

gr = 0.0
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
plot_output_flag = 1
save_output_flag = 1
plot_show_flag = 0

#search if another run with the same base name already exists
i = 0

condition = True

base_name = run_name

while condition:

    run_name = base_name + '_{0:03}'.format(i)

    backup_advanced_file = run_name + '_advanced_inp.bak'
    backup_file = run_name + '_inp.bak'

    condition = os.path.isfile(backup_file)

    i = i + 1

# create a backup file of the input parameters
shutil.copy2('run_code_simple.py', backup_file)


( h_new , h_diff , time_iter , max_angle , mean_angle  ) =       \
    ErosioNL(X,Y,h,final_time,delta_t_max,delta_t0,              \
    cr_angle,enne,k,lambda_wb,k_wb,max_nlc,max_inner_iter, res,  \
    bc, grow_rates,run_name , n_output , A_c , verbose_level,    \
    save_output_flag,plot_output_flag,plot_show_flag)
