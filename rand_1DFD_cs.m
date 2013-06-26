function [ X ] = rand_1DFD_cs( d,l )
%rand_cs_1DFD Random 1D Finite Difference Cosparse vector
%   rand_1DFD_cs(d,l) roduces a smoother vector
%   than csa_projection_1DFD(rand(d,1),l)

    L = sort(ceil(d*rand(l,1)));

    X = randn*ones(d,1);
    
    for i = 1:l-1
        X(L(i):L(i+1)-1) = randn;
    end
    
    if l > 0
        X(1:L(1)) = randn;
    end

end

