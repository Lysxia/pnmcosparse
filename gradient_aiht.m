function handle = gradient_aiht( M,y )
%gradient_aiht Gradient function
%  gradient_aiht(M,x) returns a function
%  f(x) = M'*(y-M*x)

    Mt = M';
    
    handle = @(x) Mt*(y-M*x);

end

