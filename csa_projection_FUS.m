function [ distL2,opt_pred,X,L,D ] = csa_projection_FUS( Z,l )
%csa_projection_FUS Optimal 1D Fused Lasso operator projection
%   Optimal segmentation (opt. seg.) of the signal vector Z (length d)
%   with k1 segments (k1-1 change points) and k2 non-zero entries
%   such that k1-1+k2 == k
%   k = p-l = 2d-1-l
%   (in theory, it is ... <= k, but almost never happens in practice
%   so we stick with == in the comments)
%   Precise usage ::: TODO

    d = length(Z);
    p = 2*d-1;
    k = p-l;
    
    % implicit segment indexing : index (i,j) segment is [i,j[
    % leftmost index included and rightmost index excluded
    
    X = zeros(d,1);
    
    % Trivial/nonsense cases, we don't really care here
    if d == 1
        X = Z;
        L = [];
        D = 0;
        return;
    elseif k <= 1
        X(:) = 0;
        L = [];
        D = sum(Z.^2);
        return;
    end
    
    L = zeros(k,1); % Omega*X support
    
    s = zeros(d,1);     % i-th step, s(j) = sum(Z(j'),j'=j..i)
    sumsq = zeros(d,1); % sum of the squares
    
    %   Optimal preceding index
    % i = opt_pred(j+1,c,h,omega)
    % when [i,j] is the last seg. in the opt. seg. of [1,j]
    % in k1 segments and k2 non-zero coefficients s.t. k1-1+k2 == c
    % and where the last segment [i,j] is set to zero when omega == 2
    opt_pred = ones(d,k-1,k,2);
    
    %  Optimal value on the corresponding segment [i,j]
    % Covers cases omega == 1 and omega == 2 at the same time
    opt_val = zeros(d,k-1,k);
    
    %  Minimal (L^2) distance between Z and the opt. seg. of [1,i]
    % distL2(i+1,c,k1,omega) = I_c(i,omega,k1-1) in the paper
    distL2 = zeros(d,k-1,k,2);
    
    s(1) = Z(1);
    sumsq(1) = Z(1)^2;
    
    for i = 1:d-1
        
        % omega == 1
        distL2(i+1,1:i,1,1) = sumsq(i);
        distL2(i+1,i+1:k-1,1,1) = sumsq(i) - (s(i)^2)/i;
        opt_val(i+1,i+1:k-1,1) = s(i)/i;
        
        % omega == 2
        distL2(i+1,:,1,2) = sumsq(i);
        
        var_i = sumsq(1:i) - (s(1:i).^2)./(i:-1:1)';

        [ distL2(i+1,2:k-1,2:k,2),opt_pred(i+1,2:k-1,2:k,2) ] = ...
           min(distL2(1:i,1:k-2,1:k-1,1)+repmat(var_i,[1,k-2,k-1]),[],1);
        
        s(1:i+1) = s(1:i+1) + Z(i+1);
        sumsq(1:i+1) = sumsq(1:1) + Z(i+1)^2;
    end

end

