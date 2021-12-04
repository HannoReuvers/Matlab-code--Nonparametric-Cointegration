clear variables; clc; close all;
addpath("./functions")

%% Description: Using the parametric specification test for a cubic model
% Note: We omit simulated draws for which the nonparametric estimate is
% not available for each {x= -1+0.02k; k=0,...,100}.

%% Simulation Settings
n = 5E2;
sigma = 0.2;
theta = 0.5;
KernelName = 'Epanechnikov';
h = n^(-1/3);
run('./functions/KernelDeclarations.m')
Nsim = 1E3;
MaxTrials = 10*Nsim;
PiFunction = @(x) 1;
alpha = 0.05;
rng(12345)

%% Initialize
xlist = (-1:0.02:1)';
SuccesCounter = 0;
pValue = NaN(Nsim ,1);
pValueInfeas = NaN(Nsim, 1);
tau0Estlist = NaN(Nsim, 1);
LRVEst = NaN(Nsim, 1);

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

    % Parametric estimation
    X = xt.^3;
    thetahat = (X'*yt)/(X'*X);
    resids = yt-X*thetahat;
    Eu02Est = EstResidVariance(0, xt, resids, h, MyKernel);

    % Specification test
    TnTest = TnStat(xt, resids, h, MyKernel, PiFunction, xlist);
    if TnTest>0 && isfinite(Eu02Est)
        
        % Update counter to store succesful iteration
        SuccesCounter = SuccesCounter+1;
        
        % Feasible result
        phiEst = AndrewsLongRunVariance(diff(xt));
        tau0Est = Eu02Est*muK2*(max(xlist)-min(xlist));
        TestStat = (phiEst/(sqrt(n)*h*tau0Est))*TnTest;
    
        % Infeasible result
        phi = 1;
        tau0Infeas= sigma^2*muK2*(max(xlist)-min(xlist));
        TestStatInfeas = (phi/(sqrt(n)*h*tau0Infeas))*TnTest;

        % Store p-values
        pValue(SuccesCounter) = 2*(1-normcdf(TestStat));
        pValueInfeas(SuccesCounter) = 2*(1-normcdf(TestStatInfeas));

        % Store variance estimation results
        LRVEst(SuccesCounter) = phiEst;
        tau0Estlist(SuccesCounter) = tau0Est;
    end

    % Required number of simulations ready?
    if SuccesCounter==Nsim
        break
    end
end

fprintf('\n')
disp('===== Results =====')
fprintf('Significance level (%%): %5.2f\n\n', 100*alpha);
fprintf('Empirical size [with nuisance parameter estimation] (%%): %5.2f \n', 100*mean(pValue<alpha));
fprintf('Empirical size [with known nuisance parameters] (%%): %5.2f \n\n', 100*mean(pValueInfeas<alpha));

function [fgrid] = fB(xgrid)
    % Function f_B(x)
    fgrid = xgrid.^3;
end

function [ErrorVarEst] = EstResidVariance(xpoint, x, NLSresiduals,h , K)
    ErrorVarEst = (NLSresiduals.^2'*K( (x-xpoint)/h ))/(sum( K( (x-xpoint)/h ) ));
end

function [LRV] = AndrewsLongRunVariance(y)
        % Dimensions
        T = length(y);

        % Demean input series
        y = y-mean(y);

        % <<<<< BANDWIDTH SELECTION >>>>>
        % AR(1) estimation
        yt = y(2:T);
        ylag = y(1:T-1);
        rhohat = (yt'*ylag)/(ylag'*ylag);
        residual = yt - rhohat*ylag;
        sigma2hat = var(residual);
        % alphahat1 (see Andrews (1991), Eq. (6.4))
        alphahat1num = (4*rhohat^2*sigma2hat^2)/( (1-rhohat)^6*(1+rhohat)^2 );
        alphahat1den = sigma2hat^2/(1-rhohat)^4;
        alpha1hat = alphahat1num/alphahat1den;
        % Barlett optimal bandwidth (see Andrews (1991), Eq. (5.3))
        SelectedBandwidth = (1.1447)*(alpha1hat*T)^(1/3);      
        
        % <<<<< LONG-RUN VARIANCE COMPUTATION >>>>>
        RescaledLags = (0:T-1)/SelectedBandwidth;
        BartlettWeights = (1-RescaledLags).*(RescaledLags<=1);
        LRV = y'*y/T;
        for iter = 1:(T-1)
            gammahat = y(1:(T-iter))'*y((1+iter):T)/T;
            LRV = LRV + 2*BartlettWeights(iter+1)*gammahat;
        end
end
