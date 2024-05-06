function X = rk4t(v,X,h,n,tstart)

% RK4   Runge-Kutta scheme of order 4 
%   performs n steps of the scheme for the vector field v
%   using stepsize h on each row of the matrix X
%   v maps an (m x d)-matrix to an (m x d)-matrix 

for t=tstart:h:tstart+(n-1)*h,
    k1 = v(X,t); 
    k2 = v(X + h/2*k1,t+h/2); 
    k3 = v(X + h/2*k2,t+h/2);
    k4 = v(X + h*k3,t+h);
    X = X + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end