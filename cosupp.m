function [ S ] = cosupp( x, l )
%cosupp Cosupport indices
%   cosupp(x, l) returns the indices of the l smallest elements in x
%     (in absolute value)
%     In case of equal matches, an arbitrary choice is made.

    [~,I] = sort(x);
    S = I(1:l);
end

