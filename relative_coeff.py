def relative_coeff(t,lambda_wb,k_wb):

    # RELATIVE_COEFF This function compute a non-dimensional relative
    # coefficent (0<rel_coeff<=1) multiplying the diffusion coefficient
    # k. It is based on a cumulative Weibull distribution with parameters
    # lambda_wb and k_wb

    import numpy as np

    if (lambda_wb == 0.0 ):
    
        k_rel = 1.0
        
    else:

        k_rel= np.maximum(0.001,1.0-np.exp(-(t/lambda_wb)**k_wb))

    return k_rel
