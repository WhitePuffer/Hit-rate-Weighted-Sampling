function [ keepInx, ptsNorm, th, gmm ] = GMMremove(W,GMMthreshold)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Fit two-component GMM.
ptsNorm=W;
% fprintf('Fitting dual-mode GMM for outlier removal...');

[ priors ,mu ,var ] = gmmWrapper(ptsNorm,[0.5 0.5],[min(ptsNorm) max(ptsNorm)],[max(ptsNorm)/10 max(ptsNorm)/10]);


% Compute threshold.
sig = sqrt(var);
mu = [ mu mean(ptsNorm) ];
sig = [ sig std(ptsNorm) ];
priors = [ priors 0 ];
    
    % Use GMM fitting result.
    th = 0.5*(mu(1)+mu(2));
    if priors(1)>0.8
        th = mu(1) + (priors(2)/0.2)*(th-mu(1));
    end   
% Remove low-norm points.
keepInx = ptsNorm>=GMMthreshold*th;
gmm.priors = priors;
gmm.mu = mu;
gmm.sig = sig;
end

