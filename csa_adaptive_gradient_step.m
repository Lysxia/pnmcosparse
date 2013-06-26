function [ handle ] = csa_adaptive_gradient_step( y,M )
%csa_adaptive_gradient_step Gradient Step handle for CoSparse Analysis
%   csa_adaptive_gradient_step(y,M) returns a function g(X) = Xg
%   Where Xg is obtained by gradient descent from X to lower ||y-M*Xg||

    Mt = M';
    d = size(M,2);

    function [Xg,incr] = ags( X )
        y_residual = y-M*X;
        grad = Mt*y_residual;
        a = M*grad;

        step = y_residual'*a/(a'*a);
        
        if isnan(step)
            incr(1:d) = 0;
            Xg = X;
        else
            incr = step*grad;
            Xg = X + step*grad;
        end
    end

    handle = @ags;
end