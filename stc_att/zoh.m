function [ EH,BE,ES,ZE ] = zoh( P,sol )
% ----------------------------------------------------------------------- %
%ZOH   zeroth order hold discretization
%
% ZOH(ctrl,sol) Computes the matrices (A,B,S,Z) to numerically discretize 
% a problem instance of the PDG problem about the most recent trajectory 
% stored in sol.X and sol.U. Assumes a piecewise constant thrust command.
%
% These give the first order Taylor series of the dynamics
% at the time nodes dictated by ctrl.N and correspond to a discrete time
% LTV system approximating the nonlinear dynamics.
%
% Outputs: 
%   - EH (N*Nx by N*Nx)
%   - BE (N*Nu by N*Nu)
%   - ES (N*Nx by 1)
%   - ZE (N*Nx by 1)
%
% T. Reynolds -- RAIN Lab
% ----------------------------------------------------------------------- %

% Problem dimensions
szeA    = P.Nx;
szeB1   = P.Nx;
szeB2   = P.Nu;
szeS    = P.Nx;
szeZ    = P.Nx;

% Relevant info from structs
iter    = sol.iter;
X       = sol.X(:,iter);
U       = sol.U(:,iter);
s       = sol.s(iter);
time    = P.tau;
len     = length(time)-1;
x0      = X(1:szeA);

% Set sizes for each discrete matrix
EH                  = zeros(szeA*(len+1),szeA*(len+1));
EH(1:szeA,1:szeA)   = eye(szeA);
BE                  = zeros(szeA*(len+1),szeB2*(len+1));
ES                  = zeros(szeA*(len+1),1);
ZE                  = zeros(szeA*(len+1),1);

% Initial conditions
A0  = eye(szeA);
B0  = zeros(szeB1*szeB2,1);
S0  = zeros(szeS,1);
Z0  = zeros(szeZ,1);
P0  = [x0(:); A0(:); B0(:); S0(:); Z0(:)];

for i = 1:len   
    u   = U(szeB2*(i)+1:szeB2*(i+1));
    u   = reshape(u,3,1);
    tspan = linspace(time(i), time(i+1), P.Nsub);

    F   = rk4(@(t,X)deriv(t,X,u,s,P),tspan,P0(:));

    xF  = X(szeA*i+1:szeA*(i+1)); % Restart each window with previous soln
    AF  = F(end,szeA+1:szeA*szeA+szeA);
    BF  = F(end,szeA*(szeA+1)+1:szeA*(szeA+1)+szeB1*szeB2);
    SF  = F(end,szeA*(szeA+1)+szeB1*szeB2+1:szeA*(szeA+1)+szeB1*szeB2+szeS);
    ZF  = F(end,szeA*(szeA+1)+szeB1*szeB2+szeS+1:end);
    
    % Reshape to matrices
    Ad  = reshape(AF,szeA,szeA);
    Bd  = Ad*reshape(BF,szeB1,szeB2);
    Sd  = Ad*reshape(SF,szeS,1);
    Zd  = Ad*reshape(ZF,szeZ,1);
    
    % Redefine initial condition
    P0 = [xF(:); A0(:); B0(:); S0(:); Z0(:)];
    
    % Fill up matrices
    EH(szeA*(i)+1:szeA*(i+1),szeA*(i-1)+1:szeA*i)       = Ad;
    BE(szeA*(i)+1:szeA*(i+1),szeB2*(i)+1:szeB2*(i+1))   = Bd;
    ES(szeA*(i)+1:szeA*(i+1))                           = Sd;
    ZE(szeA*(i)+1:szeA*(i+1))                           = Zd;
end


end

function DX = deriv(t,X,u,s,P)
    Nx  = P.Nx;
    x   = X(1:Nx);
    
    [f,A,B,Z] = get_vals(t,x,u,P);
    
    PHI = reshape(X(Nx+1:Nx+Nx*Nx),size(A));
       
    xdot    = s*f;
    PHI_dot = (s*A)*PHI;
    Bd_dot  = PHI\(s*B);
    Sd_dot  = PHI\f;
    Zd_dot  = PHI\(s*Z);

    DX  = [xdot(:); PHI_dot(:); Bd_dot(:); Sd_dot(:); Zd_dot(:)];
end

% ----- f expressions ----- %

function [f,A,B,Z] = get_vals(t,x,u,P)
    
    % compute jacobians    
    [A,B] = Df_att(P,x,u);
    
    % compute f
    f       = f_att(P,t,x,u);  
    
    % compute z
    Z       = - A*x - B*u;
                           
end

% ------------------------- %
