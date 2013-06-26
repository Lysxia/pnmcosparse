function [ X ] = aiht( X,params,iter,th )
%aiht Analysis Iterative Hard Thresholding
%   (This is only the skeleton of the algorithm
%   so using appropriate parameters, AHTP can also be implemented)
%   Usage ::: TODO

    if nargin < 4
        th = 0.001;
        if nargin < 3
            iter = 5;
            if nargin < 2
                throw(MException('InvalidArg:argn','Not enough arguments'));
            end
        end
    end

    % Making sure these function handles and variables are defined
    grad = params.Hgradient_step;
    cosupp = params.Hcosupp;
    projector = params.Hprojector;
    l = params.l;
    
    inner_iter = 5;
    
    Dnew = Inf;
    
    for i = 1:ceil(iter/inner_iter)
        for t = 1:inner_iter
            X_grad = grad(X);
            L = cosupp(X_grad,l);
            X = projector(X_grad,L);
        end
        Dold = Dnew;
        Dnew = norm(X-X_grad);
        if abs(Dold-Dnew) < th
            break
        end
    end

end
