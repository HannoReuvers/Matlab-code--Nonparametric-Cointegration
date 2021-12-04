function [fhat] = KernelEst(xpoint, x, y, h, K)
%% DESCRIPTION: Kernel estimator, e.g. (2.2) in Wang and Phillips (2009)
%---INPUT VARIABLE(S)---
%   (1) xpoint: point at which f(.) is estimated
%   (2) x: nonstationary regressor
%   (3) y: dependent variable
%   (4) h: bandwidth
%   (5) K: function handle for kernel
%---OUTPUT VARIABLE(S)---
%   (1) fhat: kernel estimate
    
    % Compute kernel estimate
    fhat = (y'*K( (x-xpoint)/h )) / (sum( K( (x-xpoint)/h ) ) );
    
end

