function [X,L1,L2] = csa_projection_FUS( Z,w1,w2,l )
%csa_projection_FUS Optimal sparse segmentation
%
%  Input :
%    - Signal vector Z (length d)
%    - Parameters w1,w2,l (integers)
%
%  Output :
%    - Opt. seg. of Z (see below)
%    - Support of OmegaFD * X
%      L1 : change points indices,
%      L2 : non-zero indices (as a (2,k2) matrix representing ranges of
%        non-z. segments
%
%  Returns the optimal sparse segmentation of Z with parameters w1,w2 and
%  cosparsity l.
%
%  The returned vector X is piecewise constant with :
%    - k1 change points, n_1, ..., n_k1 :
%        0 = n_0 < n_1 < ... < n_k1 < n_(k1+1) = d
%        such that for i in 0:k1 X(n_i+1:n_(i+1)) is constant
%    - k2 non-zero coefficients
%
%  subject to k1 * w1 + k2 * w2 <= p - l
%  p = w1*(d-1) + w2*d
%
%  The algorithm that is used is described (for w1=w2=1)
%  in the paper Greedy-Like Algorithms for the Cosparse Analysis Model
%  pages 9,10 and relies on dynamic programming.
%
%  Code has not been optimized for speed (yet)

  d = length(Z);
  p = w1*(d-1) + w2*d;
  k = p-l;
  
  km = k - min(w1,w2);
  
  s = zeros(d,1);   % s(n) = sum(z(1:n)) but computed efficiently
  ssq = zeros(d,1); % ...  = sum(z(1:n).^2)
  
  s(1) = Z(1);
  ssq(1) = Z(1)^1;
  
  % Data corresponding to I(n,K,1) are at opt_*(n,K)
  opt_dist = inf(d-1,km);
  opt_pred = zeros(d-1,km);
  opt_val = zeros(d-1,km);
  
  function [ dist,pred,val ] = I(L,K)
      % Last segment is zero (case omega = 0)
      dist = ssq(L);
      pred = 0;
      
      if K > w1
        for n0 = 1:L-1
            cur_dist = opt_dist(n0,K-w1) + ssq(L) - ssq(n0);
            if cur_dist < dist
                dist = cur_dist;
                pred = n0;
            end
        end
      end
      
      % Last segment is non-zero (case omega = 1)
      val = 0;
      
      nmin = ceil(L-(K-w1)/w2);
      
      % In the case of 1 segment, the associated cost would be w2*L
      % k1 = 0 and the w1 term is zero
      if (w2*L <= K)
          sn = s(L)/L;
          cur_dist = ssq(L) - (sn^2)*L;
          if cur_dist < dist
              dist = cur_dist;
              pred = 0;
              val = sn;
          end
      end
      
      % Otherwise the current segment which stops at n1 costs w1+w2*(n1-L)
      
      % This first particular case MAY reach opt_pred(~,0) (out of bounds)
      % Hence we take it out of the loop
      if (0 < nmin && nmin < L)
          sn = s(L) - s(nmin);
          cur_dist = ssq(L) - (sn^2)/(L-nmin);
          if cur_dist < dist
              dist = cur_dist;
              pred = nmin+1;
              val = sn/(L-nmin);
          end
      end
      
      for n1 = max(1,nmin+1):L-1
          sn = s(L)-s(n1);
          cur_dist = opt_dist(n1,K-w1-(L-n1)*w2) + ssq(L) - ssq(n1) - (sn^2)/(L-n1);
          if cur_dist < dist
              dist = cur_dist;
              pred = n1;
              val = sn/(L-n1);
          end
      end
  end
  
  for n = 1:d-1
      for K = 1:km
        [ opt_dist(n,K),opt_pred(n,K),opt_val(n,K) ] = ...
            I(n,K);
      end
      
      s(n+1) = s(n) + Z(n+1);
      ssq(n+1) = ssq(n) + Z(n+1)^2;
  end
  %{
  opt_dist
  opt_pred
  opt_val
  %}
  % Construct back
  
  L1 = zeros(d,1); % Change-points
  L2 = zeros(2,d); % Non-zero coefficients (in order)
  X = zeros(d,1);
  
  [ ~,L1(1),v ] = I(d,k);
  
  l1 = 1;
  l2 = 1;
  k0 = k-w1;
  
  X(L1(1)+1:d) = v;

  if v ~= 0
      L2(:,1) = [ d;L1(1)+1 ];
      l2 = l2+1;
      k0 = k0 - w2*(d+1-L1(1));
  end
  
  while L1(l1) ~= 0
    if k0 <= 0
        l1 = l1+1;
        L1(l1) = 1;
        break
    end
    L1(l1+1) = opt_pred(L1(l1),k0);
    
    v = opt_val(L1(l1),k0);
    
    X(L1(l1+1)+1:L1(l1)) = v;
    
    if v ~= 0
        L2(:,l2) = [ L1(l1);L1(l1+1)+1 ];
        l2 = l2+1;
        k0 = k0 - w1 - w2*(L1(l1) - L1(l1+1));
    else
        k0 = k0 - w1;
    end
    
    l1 = l1+1;
  end
  
  L1 = L1(l1-1:-1:1);
  L2 = reshape(L2(2*l2-2:-1:1),2,l2-1);
end