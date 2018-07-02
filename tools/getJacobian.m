function [ Df ] = getJacobian( f,x,dt )
% This function approximates the derivative of the (anonymous) function 'f'
% at the point 'x'. Tolerance 'dt' is optional.

if nargin < 3
    dt = 1e-4;
end

n = length(x);
m = length(f(x));
Df = zeros(m,n);
E = eye(m);
Ex = eye(n);

for i = 1:m
    for j = 1:n
        Df(i,j) = (1/(2*dt))*(f(x + dt.*Ex(:,j))'*E(:,i) - f(x - dt.*Ex(:,j))'*E(:,i));
    end
end

end