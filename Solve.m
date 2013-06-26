%% Solve

figure(1);
plot(1:d,x,'k',1:d,x0,'r');

figure(2);
plot(abs(x0-x));
axis([1 d 0 2]);

drawnow;

disp(struct('y_diff',D,'x_diff',norm(x0-x)));

if strcmp(strat,'asp') || strcmp(strat,'acosamp')
    x2 = asp(x0,params_FD,iter,threshold);
    x3 = asp(x0,params_generic,iter,threshold);

    %x2 = csa_projection_1DFD( x1, l_target);
    %plot(1:d,x,'k',1:d,x1,'r',1:d,x2,'b');
else
    if strcmp(strat,'aiht')
        iter = 200;
    end
    x2 = aiht(x0,params_FD,iter,threshold);
    x3 = aiht(x0,params_generic,iter,threshold);
end

D2 = norm(y-M*x2);
D3 = norm(y-M*x3);

figure(1);
plot(1:d,x,'k',1:d,x2,'r',1:d,x3,'b');

figure(2);
plot(1:d,abs(x2-x),'r',1:d,abs(x3-x),'b');
axis([1 d 0 2]);

disp(struct('y_diff_opt',D2,'x_diff_opt',norm(x2-x)));
disp(struct('y_diff_gen',D3,'x_diff_gen',norm(x3-x)));

drawnow;