function dx = f_att(P,t,x,u,ut)

q   = x(1:4);
q   = q./norm(q);
w   = x(5:7);

if( nargin == 4 )
    uu  = u;
else
    uu  = interp1(ut,u',t,P.method)';
end

dq    = 0.5*Q_mult(q,[w; 0]);
dw    = P.J\(uu - skew(w)*P.J*w);

% Costate terms %
lq  = x(8:11);
lw  = x(12:14);

dlq     = reshape(-0.5*lq'*Q_skew_star([w;0]),4,1);
temp    = -0.5*lq'*Q_skew(q);
dlw     = reshape(temp(1:3) + lw'*(P.J\(skew(w)*P.J - skew(P.J*w))),3,1);
% ------------- %

dx = [ dq; dw; dlq; dlw ];

end