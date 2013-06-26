function handle = csa_printer(x,display)

    if nargin == 1
        display = true;
    end

    d = length(x);

    function [] = output( X,D )
        figure(1);
        plot(1:d,x,'k',1:d,X,'r');
        
        figure(2);
        plot(abs(x-X));
        axis([1 d 0 1]);
        
        drawnow;
        
        if display
            disp(struct('y_diff',D,'x_diff',norm(x-X)));
        end
    end

    handle = @output;
end
