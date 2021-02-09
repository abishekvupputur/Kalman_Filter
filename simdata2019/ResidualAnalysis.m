function [avg, variance_error, Rsq] = ResidualAnalysis(Param, Input, Expected, ignore_data)
    ignore_data = 0;
    Y_estd = (Param'*Input)';
    Y = Expected(1+ignore_data:end);
    Residual = Y_estd - Y;
    avg = mean(Residual);
    variance_error = var(Residual);
    variance_signal = var(Y);
    R = corrcoef(Y_estd,Y);
    Rsq = 1 - (variance_error/variance_signal); 
end