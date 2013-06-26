function [ x1,x,y,M,Omega ] = random_inverse( d,m,p,l,noise )
    x = randn(d,1);
    M = rand_sampling_matrix(m,d);
    e = noise*(2*rand(m,1) - 1);

    % generate random tight frame with equal column norms
    if p == d
        T = randn(d);
        [Omega, ~] = qr(T);
    else
        Omega = randn(p, d);
        T = zeros(p, d);
        tol = 1e-8;
        max_j = 200;
        j = 1;
        while sum(sum(abs(T-Omega))) > tol*(p*d) && j < max_j
            j = j + 1;
            T = Omega;
            [U, ~, V] = svd(Omega);
            Omega = U * [eye(d); zeros(p-d,d)] * V';
            Omega = diag(1./sqrt(diag(Omega*Omega')))*Omega;
        end
        %disp(j);
    end

    L = ceil(p*rand(p-l,1));
    
    OmegaL = Omega(setdiff(1:p,L),:);
    
    QL = eye(d) - pinv(OmegaL) * OmegaL;
    
    x = QL * x;
    disp(norm(x));
    
    y = M*x + e;
    
    params = params_analysis(d,y,M,Omega,l,1,1,'projector');
    
    x1 = aiht(zeros(d,1),params,200,0.001);
    
    figure(1);
    plot(1:d,x,'k',1:d,x1,'r');
    
    figure(2);
    plot(abs(x-x1));
    axis([1 d 0 1]);
    
    disp(norm(x-x1));
end

