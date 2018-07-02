function [ xdot ] = Q_ode( t,x,u,ut,method )
% T. Reynolds, RAIN Lab. Updated: 8.8.17

% Quaternion ODE formula used for propagating attitude motion. Equations
% are:
%   \dot{q} = 0.5 * q \otimes w
%   \dot{w} = J^{-1} ( u - w x Jw )

q   = x(1:4);
q   = q./norm(q);
w   = x(5:7);

if nargin == 3
    uu  = u;
else
    uu  = interp1(ut,u',t,method)';
end

qdot    = 0.5*Q_mult(q,[w; 0]);
wdot    = J\(uu - skew(w)*J*w);

xdot = [qdot; wdot];

end

