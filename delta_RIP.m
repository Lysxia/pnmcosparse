function delta = delta_RIP( M,Omega,l,iter )

    [p,d] = size(Omega);
    delta = 0;
    
    for i = 1:iter
        L = setdiff(1:d-1,ceil((d-1)*rand(p-l,1)));
        OmegaL = Omega(L,:);
        Q = eye(d) - pinv(OmegaL)*OmegaL;
        N = Q*(eye(d)-M'*M)*Q;
        delta = max(delta,norm(N));
    end

end

