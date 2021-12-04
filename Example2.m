clear variables; clc; close all;
addpath("./functions")
rng(12345)

%% Description: Replication of a part of Table 1 in Wang and Phillips (2009)
% We focus on Model (B) and the estimator \hat{f}(x). The latter is
% probably preferred because the infeasible estimator is not available in
% practice and the augmented method does not always lead to improvements.
% Note: We omit simulated draws for which the nonparametric estimate is
% not available for each {x= -1+0.02k; k=0,...,100}. This requires some
% computations when the bandwidth is small (e.g. when h = n^(10/18)).

%% Simulations Settings
n = 500;
thetalist = [2; 0.2];
hpower = [-10/18; -1/2; -1/3; -1/5];
OverviewMatrix = NaN(length(thetalist)*length(hpower), 3);
for thetaiter = 1:length(thetalist)
    theta = thetalist(thetaiter);
    for hiter = 1:length(hpower)
        h = n^hpower(hiter);
        Output = SimulationRun(n, theta, h);
        [Bias, StandError, MSE] = SimulationSummary(Output);
        OverviewMatrix(4*(thetaiter-1)+hiter, :) = [Bias, StandError, MSE];
    end
end

% Lists to print nicely to screen
thetastring = {'2','0.2'};
hstring = {'n^(-10/18)','n^(-1/2)','n^(-1/3)','n^(-1/5)'};

% Print output to clean screen
clc;
fprintf('\n\n===== MODEL (B) =====\n\n')
fprintf('%15s %15s %15s %15s %15s\n', 'theta', 'h', 'Bias', 'Std', 'MSE')
for thetaiter = 1:2
    for hiter = 1:4
        fprintf('%15s %15s %15.3f %15.3f %15.3f\n', thetastring{thetaiter} ,hstring{hiter}, OverviewMatrix(4*(thetaiter-1)+hiter, 1),...
            OverviewMatrix(4*(thetaiter-1)+hiter, 2), OverviewMatrix(4*(thetaiter-1)+hiter, 3));
    end
end

% FUNCTIONS
function [fhatresults] = SimulationRun(n, theta, h)

    % Simulation parameters
    sigma = 0.2;
    Nsim = 1E3;
    MaxTrials = 1E3*Nsim;

    % Load kernel information
    KernelName = 'Epanechnikov';
    run('./functions/KernelDeclarations.m')
    

    % Simulation settings
    xlist = (-1:0.02:1)';

    % Initialize
    fhatresults = NaN(length(xlist), Nsim);
    SuccesCounter = 0;
    
    % Start simulations
    fprintf('\nStarting new simulations\n');
    for simiter = 1:MaxTrials
        % Report progress
        if mod(simiter, 5E3) == 0
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
end

function [fgrid] = fB(xgrid)
    % Function f_B(x)
    fgrid = xgrid.^3;
end

function [bias, stdev, MSE] = SimulationSummary(ResultMatrix)   
    % Dimensions
    NumberSimulations = size(ResultMatrix, 2);
    % Design info
    xlist = (-1:0.02:1)';
    ftrue = fB(xlist);
    % Simulation summary statistics
    bias = mean(ResultMatrix-ftrue*ones(1, NumberSimulations), 'all');
    MSE = mean( (ResultMatrix-ftrue*ones(1, NumberSimulations)).^2, 'all');
    stdlist = NaN(length(xlist), 1);
    for iter = 1:length(stdlist)
        stdlist(iter) = std(ResultMatrix(iter, :));
    end
    stdev = mean(stdlist);
end