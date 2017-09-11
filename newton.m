


function [x0,err,it,y,H] = newton(f,df,x0,dx,dy,maxit)
H = [x0 0 0];   % set H value
for it = 1:maxit        % do for loop
    x1 = x0-feval(f,x0)/feval(df,x0);   % set x1 value 
    err = abs(x1-x0);       % estimate the error with difference
    rerr = 2*err/(abs(x1)+dx);  % set another error
    x0 = x1;    % set x0 =x1
    y = feval(f,x0);    
    H = [H; x0 err y];                % history of iteration
    if (err<dx) || (rerr<dx) || (abs(y)<dy), break, end   % stop if the condition satisfied
end

end
