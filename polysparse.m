function [ X,L,K ] = polysparse( Z,deg,w1,w2,k )
%polysparse Optimal sparse piecewise polynomial approximation
%  Union of Subspaces projection
%  polysparse(Z,deg,w1,w2,k) returns the sparse piecewise
%  polynomial signal X which is closest to Z, with parameters (w1,w2,k).
%  (more details below)
%
%  Input :
%    - Z   : a 1D signal as a length d vector
%    - deg : an integer, max degree of the polynomial pieces that make X
%    - w1,w2,k : integer parameters
%
%  In more formal terms, the output signal X has
%    - k1 change points, 0 = n_0 < n_1 < ... < n_k1 < n_(k1+1) = d
%      such that for all l in [1,k1+1], X(n_(l-1)+1:n_l) contains
%      consecutive values of some polynomial function with degree <= deg
%    - k2 non-zero values
%      (i.e. some segments X(n_l+1:n_(l+1)) are chosen to correspond
%      in fact to the zero polynomial)
%  such that w1 * k1 + w2 * k2 <= k.
%
%  The signal X can be parameterized in terms of :
%    - The number of change points k1
%    - A set of change points      n_1 < ... < n_k1
%    - A set of polynomials        P_1 < ... < P_k1 < P_(k1+1)
%  such that X(n_(j-1)+1:n_j) = P_j(n_(j-1)+1:n_j) for j in [1,k1+1]
%
%  Denote by s_j =| 0 if P_j == 0
%                 | 1 otherwise
%
%  The signal is a solution of the following problem :
%
%    argmin[ k1 <= floor(k/w1)
%            0<n_1<...<n_k1<d
%            P_1,...,P_(k1+1) ] (
%      sum{j=1..k1+1} sum{i=n_(j-1)+1..n_j} (z(i)-P_j(i))^2
%    )
%    subject to : sum{j=1..k1+1} s_j*(n_j-n_(j-1)) <= floor((k-w1*k1)/k2)
%
%
%  Additionnal background :
%    - d-k2 is the sparsity of the signal vector X.
%      (where d is equal to length(X) and length(Z))
%      One objective of the approximation is to have as many zero values
%      as possible in the resulting signal.
%    - k is called its cosparsity, that is in this case
%      a weighted combination of sparsity and the number of change points.
%
%  Notes :
%    - it is assumed of the change points that either
%          n_(l+1) = n_l+1   (singleton segment)
%      or  n_(l+1)-n_l > deg (long enough segments
%                             to check polynomial property)
%    - the algorithm that is used here is based on dynamic programming


end