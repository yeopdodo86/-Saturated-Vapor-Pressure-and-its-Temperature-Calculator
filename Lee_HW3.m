

global ps   % we need global variable to calculate 
pslist = [1000e2;500e2;250e2];  % pressure list
for l = 1:size(pslist)  % for loop to get each pressure
    ps = pslist(l);
end
t0=60.00;   % initial
ds=1.e-5;   % delta s
dt=1.e-5;   % delta t
maxit=100;  % upper bound
[x0,err,it,y,H] = newton(@svpd,@dsdt,t0,dt,ds,maxit); % call newton to calculate
[c,yc,err,H] = bisect(@svpd,60,120,1.e-5); % call bisection to calculate

%%
%[H]=boiling_point_temp(pres15_list);

% [~,~,~,H] = bisect(f,a,b,delta);
%disp(H);

 

 
function [H] = boiling_point_temp(press_list)

[new_press] = satprs(press_list);
% new_press = press_list

a = 0.00;   % initial value for left end
b = 0.11;   % initial value for right end
  % ball function
delta = 1.e-5;  % x tolerance
[~,~,~,H] = bisect1(f,a,b,delta);
disp(H);
 
 
 
 
 
 
f = inline('x^3-0.165*x^2+3.993e-4');   % ball function
df = inline('3*x^2-0.33*x');            % function derivative
x0 = 0.05;  % initial value
dx = 1.e-5; % x tolerance
dy = 1.e-15;    % f tolerance
maxit = 100;    % max iterations
[x0,err,it,y,H] = newton1(f,df,x0,dx,dy,maxit);
fprintf('\n--it-----------x0----------err------------y\n');
for k = 2:it+1
     fprintf('%4d %12.4g %12.4g %12.4g \n',k-1,H(k,:));
end
 
 
end
 
 
 
 function [x0,err,it,y,H] = newton1(f,df,x0,dx,dy,maxit)
H = [x0 0 0];
for it = 1:maxit    
    x1 = x0-feval(f,x0)/feval(df,x0);
    err = abs(x1-x0);
    rerr = 2*err/(abs(x1)+dx);
    x0 = x1;
    y = feval(f,x0);
    H = [H; x0 err y];                % history of iteration
    if (err<dx) || (rerr<dx) || (abs(y)<dy), break, end
end

end
 
 
 
 
function [c,yc,err,H] = bisect1(f,a,b,delta)
ya = feval(f,a); yb = feval(f,b);
H = [a b ya];
if ya*yb > 0, end
maxit = 1 + round((log(b-a)-log(delta))/log(2));
for it = 1:maxit
   c  = (a+b)/2; yc = feval(f,c);
   if (yc == 0), a = c; b = c;
   elseif (yb*yc > 0), b = c; yb = yc;
   else
       a = c; 
       ya = yc;
   end
   H = [H;a b yc];
   if b-a < delta, break, end
end
c  = (a+b)/2; yc = feval(f,c);
err = abs(b-a)/2;

end

 
 
 
 
 
 
 
 
 
 
 


