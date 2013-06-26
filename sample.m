function [ D,var ] = sample( d,m,l,algo,param,noise,iter,draw )
    D = 0;
    var = 0;

    e = 0;
    t = 10*iter;

    while e < iter && t > 0
        t = t-1;
        
        x = rand_1DFD_cs(d,d-1-l);
        M = rand_sampling_matrix(m,d);
        y = M*x+noise*rand(m,1);

        x1 = zeros(d,1);
        
        par = param(y,M,l);

        x1 = algo(x1,par);

        if draw
            figure(2);
            plot(abs(x1-x),'b');
            axis([1 d 0 1]);
            drawnow;
        end

        N = norm(x-x1)/norm(x);

        if ~isnan(N)
            D = D + N;
            var = var + N^2;
            e = e + 1;
        end
    end

    D = D/iter;
    var = var/iter - D^2;

end

