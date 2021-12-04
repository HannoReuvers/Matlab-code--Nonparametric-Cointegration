clear variables; clc; close all;
addpath("./functions")

%% Description: Replication of Figure 4 in Wang and Phillips (2009)
% Note 1: The text suggests to compute the 95% "estimation band" such that
% the interval [f(x_j)-\delta, f(x_j)+\delta] contains 95% of the
% simulation draws. For computational convenience, we use quantiles.
% Note 2: We omit simulated draws for which the nonparametric estimate is
% not available for each {x= -1+0.02k; k=0,...,100}.

%% Simulation Settings
n = 500;
sigma = 0.2;
theta = 2;
KernelName = 'Epanechnikov';
h = n^(-1/3);
run('./functions/KernelDeclarations.m')
Nsim = 5E3;
MaxTrials = 10*Nsim;
rng(12345)

%% Initialize
xlist = (-1:0.02:1)';
SuccesCounter = 0;
fhatresults = NaN(length(xlist), Nsim);

%% Start simulations
fprintf('\nStarting simulations\n');
for simiter = 1:MaxTrials

    % Report progress
    if mod(simiter, 5E2) == 0
        fprintf('\tIteration %5d out of %5d: %5d out of %5d results \n', simiter, MaxTrials, SuccesCounter, Nsim);
    end

    % Generate data
    epsit = normrnd(0, 1, [n 1]);
    lambdat = normrnd(0, 1, [n 1]);
    xt = cumsum(epsit);
    ut = (lambdat+theta*epsit)/sqrt(1+theta^2);
    yt = fB(xt) + sigma*ut;

    % Kernel estimator
    fhatlist = NaN(length(xlist), 1);
    for xiter = 1:length(xlist)
        xpoint = xlist(xiter);
        fhatlist(xiter) = KernelEst(xpoint, xt, yt, h, MyKernel);
    end

    % Check for complete nonparametric estimate
    if sum(isnan(fhatlist))==0
        SuccesCounter = SuccesCounter+1;
        fhatresults(:, SuccesCounter) = fhatlist;
    end

    % Required number of simulations ready?
    if SuccesCounter==Nsim
        break
    end
end

% Compute estimation bands
[BLow, BUp] = EstimationBands(fhatresults, 0.05);

figure(1)
truef = fB(xlist);
plot(xlist, truef, 'r', 'LineWidth', 2)
hold on
plot(xlist, mean(fhatresults,2), '--k', 'LineWidth', 2)
plot(xlist, BLow, ':k', 'LineWidth', 2)
plot(xlist, BUp, ':k', 'LineWidth', 2)
hold off
axis([-1 1 -1.2 1.2])
box on

function [fgrid] = fB(xgrid)
    % Function f_B(x)
    fgrid = xgrid.^3;
end

function [BandLow, BandUp] = EstimationBands(InputMatrix, alpha)
    % Dimensions
    NumberPoints = size(InputMatrix, 1);
    % Quantiles
    BandLow = NaN(NumberPoints, 1);
    BandUp = NaN(NumberPoints, 1);
    for iter = 1:NumberPoints
        BandLow(iter) = quantile(InputMatrix(iter, :), alpha/2);
        BandUp(iter) = quantile(InputMatrix(iter, :), 1-alpha/2);
    end
end

    
