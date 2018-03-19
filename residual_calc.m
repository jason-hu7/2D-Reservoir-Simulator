function [residual] = residual_calc(T, P_n_1, R_D, R_w)
    
residual = T * P_n_1 - R_D + R_w;


end