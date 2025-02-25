function [res] = BC_jaco_MEE(CS_i,PD)

% input:  
%         CS_i, values for costates (produced by the solver)
%         PD,   problem data structure PD
% output: 
%         residual vector

% extract needed parameters from PD
options_ode = PD.options_ode;
S_i         = PD.S_i;
S_f         = PD.S_f;
tspan       = PD.tspan;

% Partial of contraints at t_f by partial of states and costates at t_f
%P1          = PD.P1;
% Partial of states and costates at t_0 by partial of costates at t_0
%P2          = PD.P2;

T  = PD.T;
c = PD.c;
eta = PD.eta;
rho = PD.rho;

% define the vector of initial conditions, which is a 14 by 1 vector
X_i = [ S_i;   % states
       CS_i;];   % costates  
       %stm_i]; % STM
   
% propagate the differential equations   
%[time,X_p] = ode45(@EOM_jaco,tspan,X_i,options_ode,PD); 
%[time,X_p] = ode45_fast_mex(tspan,X_i,c, T, eta, rho); 
[~,X_p] = fastode45_MEE_mex(tspan,X_i,c, T, rho, eta); 

% extract the final solution vector 
Z_f     = X_p; % making into a column vector

%Z_f    = X_p(:);
p_f  = Z_f(1);
f_f  = Z_f(2);
g_f  = Z_f(3);
h_f = Z_f(4);
k_f = Z_f(5);
L_f = Z_f(6);
m_f = Z_f(7);
lp_f  = Z_f(8);
lf_f  = Z_f(9);
lg_f  = Z_f(10);
lh_f = Z_f(11);
lk_f = Z_f(12);
lL_f = Z_f(13);
lm_f = Z_f(14);
%stm_f_vec1 = Z_f(15:end);

%stm_f_mat = reshape(stm_f_vec1,[14,14]);

% Final constraints
C1 = p_f - S_f(1);
C2 = f_f - S_f(2);
C3 = g_f - S_f(3);
C4 = h_f - S_f(4);
C5 = k_f - S_f(5);
C6 = L_f - S_f(6);

% Derived using transversality condition
C7 = lm_f + 1;

res = [C1;C2;C3;C4;C5;C6;C7];

%jaco = P1*stm_f_mat*P2;

end