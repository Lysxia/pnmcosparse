function [ D,var ] = csa_phasediag( d, delta_range,rho_range,algo_handle,param_handle,noise,iter )
%csa_phasediag CoSparse Analysis Phase Diagram
%   csa_phasediag(d,handle) where handle(y,M,xinit)
%   iterates starting from xinit using some specified algorithm
%   to approximate some value x such that y = Mx+e
    if nargin < 7
        throw(Mexception('InvalidArg:argn','Not enough arguments'));
    end
    
    D = zeros(rho_range,delta_range,iter);
    var = zeros(rho_range,delta_range,iter);
    
    parfor p = 1:delta_range*rho_range*iter
        rho = rem(p-1,rho_range) + 1;
        q = floor((p-1)/rho_range);
        delta = rem(q,delta_range) + 1;
            
        m = ceil(delta*d/delta_range);
        l = ceil(d - rho*m/rho_range);

        [D(p),var(p)] = sample(d,m,l,algo_handle,param_handle,noise,10,false);

        disp({delta rho D(p) m l});
    end
    
    D = mean(D,3);
    var = mean(var,3);

    D(:) = 1-min(D(:),1);

end
