function X = csa_projection_1DFD_support( X,L )
    
    l = length(L);
    d = length(X);
    
    if l == 0
        X(:) = mean(X(:));
        return
    end

    L = sort(L);

    X(1:L(1)) = mean(X(1:L(1)));
    for i = 1:l-1
        X(L(i)+1:L(i+1)) = mean(X(L(i)+1:L(i+1)));
    end
    X(L(l)+1:d) = mean(X(L(l)+1:d));

end

