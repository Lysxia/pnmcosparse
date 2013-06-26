function [ X,L ] = asp( X,params,iter,threshold )
%asp Analysis Subspace Pursuit
%   (This is only the skeleton of the algorithm
%   so using appropriate parameters, ACoSaMP can also be implemented)
%   Usage ::: TODO
    
    if nargin < 4
        threshold = 0.001;
        if nargin < 3
            iter = 5;
            if nargin < 2
                throw(MException('InvalidArg:argn','Not enough arguments'));
            end
        end
    end

    % Making sure these function handles and variables are defined
    residue = params.Hresidue;
    cosupp = params.Hcosupp;
    intersect = params.Hintersect;
    tmp_estim = params.Htmp_estim;
    projector = params.Hprojector;
    l = params.l;
    al = params.al;
    L = params.L;
    y = params.y;
    M = params.M;
    
    inner_iter = 3;
    th = inner_iter * threshold;
    
    Dnew = Inf;
    
    for i = 1:ceil(iter/inner_iter)
        for j = 1:inner_iter % Giving ourselves some more margin
            LD = cosupp(residue(X),al);
            Lt = intersect(L,LD);
            X_tmp = tmp_estim(Lt);
            L = cosupp(X_tmp,l);
            X = projector(X_tmp,L);
        end
        Dold = Dnew;
        Dnew = norm(y-M*X);
        if abs(Dold-Dnew) < th
            break
        end
    end

end

