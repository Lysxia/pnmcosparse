function [ params ] = params_1DFD( d,y,M,l,a,lambda,projector_type )
%params_1DFD Parameters for aiht and asp functions (1DFD operator)

    m = length(y);
    [mM,dM] = size(M);
    
    if m ~= mM || d ~= dM
        throw(MException('InvalidArg:sizes','Sizes do not match'));
    end

    params = struct();
    
    params.l = l;
    params.a = a;
    params.al = ceil(a*l);
    params.d = d;
    params.y = y;
    params.M = M;
    params.lambda = lambda;                 % Purely informative
    params.projector_type = projector_type; %
    params.p = d-1;                         %
    
    Mt = M';
    set1top = 1:d-1;
    Omega = fdamatrix(d);
    
    params.Hgradient_step = csa_adaptive_gradient_step(y,M);
    
    function XL = XLstruct(X,L)
        XL = struct('X',X,'L',L);
    end
    
    function XL = cosupp(X,l)
        [X,L] = csa_projection_1DFD(X,l);
        XL = XLstruct(X,L);
    end

    function XL = intersect(L0,L1)
        XL = XLstruct(L0.X,union(L0.L,L1.L));
    end

    function X = relaxed_min(L)
        X = csa_projection_1DFD_support(csa_rmin(y,M,Mt,Omega(setdiff(set1top,L.L),:),lambda),L.L);
    end

    params.Hresidue = @(X) Mt*(y-M*X);
    params.Hcosupp = @cosupp;
    params.Hintersect = @intersect;
    params.Htmp_estim = @relaxed_min;
    
    if strcmp(projector_type,'projector')
        params.Hprojector = @(~,XL) XL.X;
    elseif strcmp(projector_type,'estimate')
        params.Hprojector = @(~,XL) relaxed_min(XL);
    else
        throw(MException(...
            'InvalidArg:projector_type','Unknown projector type'));
    end

    params.L = XLstruct(zeros(d,1),[]);
    
end

