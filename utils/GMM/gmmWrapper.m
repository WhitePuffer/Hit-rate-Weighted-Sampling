function [ prior, mu, sig] = gmmWrapper(x, prior, mu, sig)
%[ prior mu sig ] = gmmWrapper(x,prior,mu,sig)
%Wrapper code for Sylvain Calinon's MATLAB GMM fitting program.
%
%x (1 x n)
%    Input data.
%prior (1 x 2)
%    Initial weight of the two means.
%mu (1 x 2)
%    Initial location of the two means.
%sig (1 x 2)
%    Initial std. dev. of the two means.
%prior (1 x 2)
%    Weights of the two means after GMM fitting.
%mu (1 x 2)
%    Location of the two means after GMM fitting.
%sig (1 x 2)
%    Std. dev. of the two means after GMM fitting.
%
%Copyright (c) 2009 Tat-Jun Chin
%School of Computer Science, The University of Adelaide, South Australia
%http://www.cs.adelaide.edu.au/~tjchin
%
%This program is part of the package that implements the paper:
%T.-J. Chin, H. Wang and D. Suter
%Robust Fitting of Multiple Structures: The Statistical Learning Approach
%In Proc. Int. Conf. on Computer Vision 2009, Kyoto, Japan
%
%The program is free for non-commercial academic use. Any commercial use
%is strictly prohibited without the author's consent. Please acknowledge
%the authors by citing the above paper in any academic publications that
%have made use of this program or part of it.

% addpath('GMM-GMR-v1.2');

var = zeros(1,1,2);
var(1,1,1) = sig(1);
var(1,1,2) = sig(2);

[ prior, mu, var] = EM(x,prior,mu,var);

sig = [ var(1,1,1) var(1,1,2) ];

end
