function [ X,L,D ] = csa_projection_1DFD( Z,l )
%csa_projection_1DFD Optimal 1D Finite Difference operator projection
%   Optimal segmentation (opt. seg.) of the signal vector Z (length d)
%   with d-l segments
%   Precise usage ::: todo

    d = length(Z);
    
    k = d-l; % max. number of segments
             % implicit segment indexing : index (i,j) segment is [i,j-1]
             % or [i,j[, leftmost index is included, rightmost is excluded
    
    X = Z(:);
    
    % Trivial/nonsense cases
    if k <= 1
        X(:) = mean(Z);
        L = [];
        D = var(Z);
        return;
    end
    
    L = zeros(k-1,1);

    s = zeros(d,1);     % during the i-th step, s(j) = sum(Z(j'),j'=j..i)
    sumsq = zeros(d,1); % sum of the squares
                        
    s(1) = Z(1);
    sumsq(1) = Z(1)^2;
    
    % Optimal preceding index, such that i = opt_pred(j+1,l) when [i,j]
    %   is the last segment in the opt. seg. of [1,j] in l segments
    opt_pred = ones(d,k-1);
    
    % Optimal value for X on the corresponding segment [i,j]
    % (same indexing convention)
    opt_val = zeros(d,k-1);
    
    % distL2(i+1,l) : minimal L^2 distance between Z and
    %   the opt. seg. of [1,i] with at most l segments
    distL2 = zeros(d,k-1);
    
    for i = 1:d-1
        % Optimal segmentation of [1,i] ([1,i+1[)
        
        distL2(i+1,1) = sumsq(1) - (s(1)^2)/i; % solution for 1 segment
        
        [ distL2(i+1,2:k-1),opt_pred(i+1,2:k-1) ] = ...
            min( distL2(1:i,1:k-2) + repmat(sumsq(1:i) - (s(1:i).^2)./(i:-1:1)', 1,k-2) );
        
        opt_val(i+1,:) = s(opt_pred(i+1,:))'./(i+1-opt_pred(i+1,:));
        
        s(1:i+1) = s(1:i+1) + Z(i+1);
        sumsq(1:i+1) = sumsq(1:i+1) + Z(i+1) ^ 2;
    end
    
    [ D,L(1) ] = min(distL2(:,k-1) + sumsq - (s.^2)./(d:-1:1)');
    
    X(L(1):d) = s(L(1))/(d+1-L(1));
    
    l = 1;
    
    while L(l) ~= 1
        L(l+1) = opt_pred(L(l),k-l);
        X(L(l+1):L(l)-1) = opt_val(L(l),k-l);
        l = l+1;
    end
    
    L = L(1:l-1) - 1;

end
