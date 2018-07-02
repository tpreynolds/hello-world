function [ A,B ] = Df_att( P,x,u )

q   = x(1:4);
w   = x(5:7);
lq  = x(8:11);
lw  = x(12:14);

% Kinematic terms
dfq_dq  = 0.5*Q_skew_star([w;0]);
dfq_dw  = 0.5*Q_skew(q);
dfq_dw  = dfq_dw(:,1:3); % delete scalar column
dfq_dlq = zeros(4,4);
dfq_dlw = zeros(4,3);

% Dynamic terms
dfw_dq  = zeros(3,4);
dfw_dw  = -P.J\(skew(w)*P.J - skew(P.J*w));
dfw_dlq = zeros(3,4);
dfw_dlw = zeros(3,3);

% Costate terms
dflq_dq     = zeros(4,4);
dflq_dw     = 0.5*Q_skew(lq);
dflq_dw     = dflq_dw(:,1:3);
dflq_dlq    = -0.5*Q_skew_star([w;0]);
dflq_dlw    = zeros(4,3);

dflw_dq     = -0.5*Q_skew(lq);
dflw_dq     = dflw_dq(1:3,:); % delete scalar column
dflw_dw     = zeros(3,3);
dflw_dlq    = -0.5*Q_skew(q);
dflw_dlq    = dflw_dlq(1:3,:);
dflw_dlw    = P.J\(skew(w)*P.J - skew(P.J*w));

% Input terms
dfq_du      = zeros(4,3);
dfw_du      = inv(P.J);
dflq_du     = zeros(4,3);
dflw_du     = zeros(3,3);

A   = [ dfq_dq  dfq_dw  dfq_dlq     dfq_dlw;
        dfw_dq  dfw_dw  dfw_dlq     dfw_dlw;
        dflq_dq dflq_dw dflq_dlq    dflq_dlw;
        dflw_dq dflw_dw dflw_dlq    dflw_dlw];
    
B   = [ dfq_du; 
        dfw_du;
        dflq_du;
        dflw_du ];    


end
