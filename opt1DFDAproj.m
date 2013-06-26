function [ X,L,D ] = opt1DFDAproj( Z,l )
%opt1DFDAproj Optimal 1D Finite Difference Projection
%   Optimal segmentation (opt. seg.) of the signal vector Z (length d)
%   with d-l segments

    Z = Z(:);

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

    s = zeros(d,1);     % at the i-th step, s(j) = sum(Z(j'),j'=j..i)
    sumsq = zeros(d,1); % sum of the squares
                        % These are updated ASAP in the loop, see below
    
    % Optimal preceding index, such that i = opt_pred(j+1,l) when [i,j]
    %   is the last segment in the opt. seg. of [1,j]
    opt_pred = zeros(d,k-1);
    
    % Optimal value for X on the corresponding segment [i,j]
    % (same indexing convention)
    opt_val = zeros(d,k-1);
    
    % distL2(i+1,l) : minimal L^2 distance between Z and
    %   the opt. seg. of [1,i] with at most l segments
    distL2 = zeros(d,k-1);
    
    for i = 1:d-1
        % Optimal segmentation of [1,i] ([1,i+1[)
        
        s(1:i) = s(1:i) + Z(i);
        sumsq(1:i) = sumsq(1:i) + Z(i) ^ 2;
        
        distL2(i+1,1) = sumsq(1)/i - (s(1)/i) ^ 2; % solution for 1 segment
        
        [ distL2(i+1,2:k-1),opt_pred(i+1,2:k-1) ] = ...
            min( distL2(1:i,1:k-2) + repmat((sumsq(1:i)./(i:-1:1)) - (s(1:i)./(i:-1:1)) .^ 2, 1,k-2) );
        
        opt_val(i+1,:) = s(opt_pred(i+1,:))./(i+1-opt_pred(i+1,:));
    end
    
    [ D,L(1) ] = min(distL2(1:d,k-1) + sumsq(1:d)./(d:-1:1) - (s(1:d)./(d:-1:1)) .^ 2);
    
    X(L(1):d) = 
    
    l = 1;
    
    while L(l) ~= 0
        L(l+1) = opt_pred(L(l),k-l);
        X(L(l+1):L(l)-1) = 

end
