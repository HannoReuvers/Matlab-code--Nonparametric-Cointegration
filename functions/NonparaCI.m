function [confLow, confUp] = NonparaCI(xpoint, x, y ,h , K, muK1, muK2, alpha)
%% DESCRIPTION: Confidence interval for nonparametric regression function, see page 1918 in Wang and Phillips (2009)
%---INPUT VARIABLE(S)---
%   (1) xpoint: point at which f(.) is estimated
%   (2) x: nonstationary regressor
%   (3) y: dependent variable
%   (4) h: bandwidth
%   (5) K: function handle for kernel
%   (6) muK1: integral over kernel function
%   (7) muK2: integral over squared kernel function
%   (8) alpha: confidence interval will have (asymptotic) coverage 1-alpha
%---OUTPUT VARIABLE(S)---
%   (1) confLow: lower bound of confidence interval
%   (2) confUp: upper bound of confidence interval

    % Compute kernel estimate
    fEst = KernelEst(xpoint, x, y, h, K);

    % Quantities to determine CI width
    temp = sum( K( (x-xpoint)/h ));
    sigma2hat = ((y-fEst).^2'*K( (x-xpoint)/h ))/(sum( K( (x-xpoint)/h ) ) );

    % Lower and upper bound of pointwise confidence interval
    confLow = fEst - norminv(1-alpha/2)*sqrt( (sigma2hat*muK2/muK1)/temp  );
    confUp = fEst + norminv(1-alpha/2)*sqrt( (sigma2hat*muK2/muK1)/temp  );
end

