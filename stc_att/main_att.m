% Main file for STC based attitude reorientation
clear;
addpath(genpath('../tools/'))

P.Nx    = 14;
P.Nu    = 3;
P.N     = 80;
P.Nsub  = 10;
P.method = 'previous';

P.q0  = Q_rand(5);
P.qf  = [ 0.0; 0.0; 0.0; 1.0 ];
P.w0  = zeros(3,1);
P.wf  = zeros(3,1);
P.x0    = [ P.q0; P.w0 ];
P.xf    = [ P.qf; P.wf ];

P.tf  = 2.0;
P.tau = linspace(0,1,P.N);

P.J = diag([1.0,1.0,2.0]); % symmetric for now

% Initial guess
qt  = Q_SLERP(P.q0,P.qf,P.tau);
wt  = reshape(repmat(P.w0,P.N,1),3,P.N);
for k = 1:P.N
    sol.X(:,k)  = [ qt(:,k); wt(:,k); zeros(7,1) ];
    sol.U(:,k)  = zeros(P.Nu,1);
end
sol.X   = vec(sol.X);
sol.U   = vec(sol.U);
sol.s   = P.tf;

P.iter_max = 10;
P.w_nu  = 1e4;
P.w_eta = 1e-1;
P.w_s   = 1e-2;

for k = 1:P.iter_max
   sol.iter = k;
   
   % Linearize
   [EH,BE,ES,ZE] = zoh(P,sol);
   
   XL = reshape(sol.X(:,k),P.Nx,P.N);
   sL = sol.s(k);
   
   % Solve
   cvx_begin quiet
    cvx_precision('low')
    
    % VARIABLES
    variables x(P.Nx*P.N) u(P.Nu*P.N)
    variable s nonnegative
    variables eta(P.N,1) eta_s v(P.Nx*P.N)
    
    % COST FUNCTION
    minimize( s + P.w_nu*norm(v,1) + P.w_eta*norm(eta) + P.w_s*norm(eta_s) )
    
    % CONSTRAINTS
    subject to
    x(1:7)   == P.x0;
    x           == EH*x + BE*u + ES*s + ZE + v;
    x(P.Nx*(P.N-1)+1:P.Nx*(P.N-1)+7) == P.xf;
    
    (s-sL)'*(s-sL) <= eta_s;
    
    for n = 1:P.N        
       lwl = XL(12:14,n); 
       un   = u(P.Nu*(n-1)+1:P.Nu*n);
       norm(un,inf) <= 1.0; 
              
%        for m = 1:3
%            -min(-lwl(m),0)*(ul(m)-1.0) == 0;
%            -min(lwl(m),0)*(ul(m)+1.0)  == 0;
%        end
    end
    
    cvx_end
    
    sol.X(:,k+1) 	= x;
    sol.U(:,k+1)    = u;
    sol.s(k+1)      = s;
    
    fprintf('Iter: %d |',k)
    fprintf(' HoG: %2.2e |',norm(v,1))
%     fprintf(' TR: %2.2e |',norm(eta))
    fprintf(' t_f: %2.2f\n',s)
    
    if( k > 1 )
        if( (norm(v,1) < 1e-8) && (abs(sol.s(k+1)-sol.s(k)) < 1e-3) )
            fprintf('Converged.\n')
            break;
        end
    end
end

%% Results
close all

s   = sol.s(:,end);
X   = reshape(sol.X(:,end),P.Nx,P.N);
U   = reshape(sol.U(:,end),P.Nu,P.N);
t   = linspace(0,s,P.N);

for k = 1:P.N
    uk = U(:,k);
    xk = X(:,k); 
    
    trj.u_mag(k) = norm(uk,inf);
end

figure(1)
for k = 1:3
subplot(3,1,k), hold on, grid on
plot(t,U(k,:),'Linewidth',1)
xlabel('Time [s]')
xlim([0 s])
end
% plot(t,trj.u_mag,'r--','LineWidth',1)
% legend('u_x','u_y','u_z','||u||')

figure(2)
subplot(2,1,1), hold on, grid on
plot(t,X(1:4,:),'LineWidth',1)
xlim([0 s])
title('Quaternion')
subplot(2,1,2), hold on, grid on
plot(t,X(5:7,:),'LineWidth',1)
title('Angular Velocity')
xlim([0 s])

figure(3)
subplot(2,1,1), hold on, grid on
plot(t,X(8:11,:),'LineWidth',1)
xlim([0 s])
title('Costate: Quaternion')
subplot(2,1,2), hold on, grid on
plot(t,X(12:14,:),'LineWidth',1)
title('Costate: Angular Velocity')
xlim([0 s])
