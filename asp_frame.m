function [ X ] = asp_frame( y,M,l,al,Xinit,p,handle,max_iteration )
%asp_frame ASP and ACoSaMP frame function
    
    X = Xinit;
    
    if isfield(handle,{'cosupp','proj1','proj2'})
        cosupp = handle.cosupp;
        proj1 = handle.proj1;
        proj2 = handle.proj2;
    else
        throw(MException('InvalidArg:MissingHandle','Handle parameter incomplete'));
    end
    
    % Constant
    Mt = M';
    
    L = 1:p;
    
    for t = 1:max_iteration
        y_resid = y-M*X;
        L1 = cosupp(Mt*y_resid,al);
        L = intersect(L,L1);
        X = proj1(L);
        L = cosupp(X,l);
        X = proj2(X,L);
    end

end

