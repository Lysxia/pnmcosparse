function [ params ] = params_1DFD_subopt( d,y,M,l,lambda,a,projector_type,reference_iter )
%params_1DFD_subopt Parameters for aiht and asp functions
%  Using suboptimal projection

    m = length(y);
    [mM,dM] = size(M);
    
    if m ~= mM || d ~= dM
        throw(MException('InvalidArg:sizes','Sizes do not match'));
    end

    params = struct();
    
    params.l = l;
    params.a = a;
    params.al = ceil(a*l);
    params.d = d;                           %
    params.y = y;                           %
    params.M = M;                           %
    params.lambda = lambda;                 % Purely informative
    params.projector_type = projector_type; %
    params.p = d-1;                         %
    
    if nargin == 8
        params.r_iter = reference_iter;
    end
    
    Mt = M';
    set1top = 1:d-1;
    Omega = fdamatrix(d);
    
    params.Hgradient_step = csa_adaptive_gradient_step(y,M);
    
    function L = cosupp(X,l)
        [ ~,L ] = sort(abs(X(1:d-1)-X(2:d)));
        L = L(l+1:d-1);
    end

    function L = intersect(L0,L1)
        L = union(L0,L1);
    end

    function X = relaxed_min(L)
        X = csa_projection_1DFD_support(csa_rmin(y,M,Mt,Omega(setdiff(set1top,L),:),lambda),L);
    end

    params.Hresidue = @(X) Mt*(y-M*X);
    params.Hcosupp = @cosupp;
    params.Hintersect = @intersect;
    params.Htmp_estim = @relaxed_min;
    
    if strcmp(projector_type,'projector')
        params.Hprojector = @csa_projection_1DFD_support;
    elseif strcmp(projector_type,'estimate')
        params.Hprojector = @(~,L) relaxed_min(L);
    else
        throw(MException(...
            'InvalidArg:projector_type','Unknown projector type'));
    end

    params.L = [];
    
end

