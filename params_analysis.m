function [ params ] = params_analysis( d,y,M,Omega,l,a,lambda,projector_type )
%params_analysis Parameters for aiht and asp functions

    m = length(y);
    [mM,dM] = size(M);
    [p,dO] = size(Omega);
    
    if m ~= mM || d ~= dM || dO ~= d
        throw(MException('InvalidArg:sizes','Sizes do not match'));
    end

    params = struct();
    
    params.l = l;
    params.al = ceil(a*l);
    params.d = d;
    params.y = y;
    params.M = M;
    params.p = p;
    params.lambda = lambda;                 % Purely informative
    params.projector_type = projector_type; %
    
    Mt = M';
    
    function L = cosupp(X,l)
        z = Omega*X;
        [ ~,L ] = sort(abs(z));
        for k = l+1:p
            if z(L(k)) ~= 0
                break
            end
        end
        L = L(1:max(k-1,l));
    end

    function X = relaxed_min(L)
        X = csa_rmin(y,M,Mt,Omega(L,:),lambda);
    end

    function X = project(X,L)
        OmegaL = Omega(L,:);
        QL = eye(d) - pinv(OmegaL)*OmegaL;
        X = QL*X;
    end

    % params.Hprojector
    if strcmp(projector_type,'projector')
        params.Hprojector = @project;
    elseif strcmp(projector_type,'estimate')
        params.Hprojector = @(~,L) relaxed_min(L);
    else
        throw(MException(...
            'InvalidArg:projector_type','Unknown projector type'));
    end
    
    params.Hgradient_step = csa_adaptive_gradient_step(y,M);
    params.Hresidue = @(X) Mt*(y-M*X);
    params.Hcosupp = @cosupp;
    params.Hintersect = @intersect;
    params.Htmp_estim = @relaxed_min;
    params.L = 1:p;

end