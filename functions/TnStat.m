function [Tn] = TnStat(x, NLSresid, h, K, Pi, testgrid)
%% DESCRIPTION: Tn statistic for specification test, see page 366 in Wang and Phillips (2016)
%---INPUT VARIABLE(S)---
%   (1) x: nonstationary regressor
%   (2) NLSresid: residuals from the parametric model estimated by NLS
%   (3) h: bandwidth
%   (4) K: function handle for kernel
%   (5) Pi: function handle for pi-function
%   (6) testgrid: grid of values to approximate integral
%---OUTPUT VARIABLE(S)---
%   (1) Tn: test statistic

    % Dimensions
    n = length(x);
    gridpoints = length(testgrid);

    % Quadratic expression in integral for Tn
    QuadraticFormGrid = NaN(length(testgrid), 1);
    for xiter = 1:gridpoints
        % Select x
        xpoint = testgrid(xiter);
        % Compute quadratic term in integral of Tn
        QuadraticFormAtx = 0;
        for k = 1:n
            QuadraticFormAtx = QuadraticFormAtx+K( (x(k)-xpoint)/h )*NLSresid(k);
        end
        QuadraticFormGrid(xiter) = QuadraticFormAtx;
    end

    % pi(x) at values in xgrid
    PiGrid = NaN(length(testgrid), 1);
    for xiter = 1:gridpoints
        PiGrid(xiter) = Pi( testgrid(xiter) );
    end

    % Tn
    Tn = sum( (QuadraticFormGrid(1:end-1).^2).*PiGrid(1:end-1).*diff(testgrid) );
end

