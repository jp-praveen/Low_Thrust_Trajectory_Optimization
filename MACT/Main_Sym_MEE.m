%{
Author: Praveen
Department of Aerospace Engineerig
Auburn University
Date modified: 8/19/2020
Output: This function generates a file for the set of state-costates that
will be used as the function hanlde for ode45 or for fastode
%}
clear
clc
NofEq = 7; % number of state differential equations


syms p f g h k L m ur ut un w s_sq Thr c mu rho delta SI2CAN t

X = [p;f;g;h;k;L;m];     
U = [ur;ut;un];  % control unit direction vector

%% Create co-state vector
P = cell(NofEq,1);
for i=1:NofEq
    P{i} = sprintf('P%d',i);
end
P = P(:); % now x is a 7-by-1 vector
P = sym(P,'real');

%% States dynamics
mu = 1;

% thurst acceleration
Ta = Thr/m*SI2CAN*delta*U;
w = 1 + f*cos(L) + g*sin(L);  
s_sq = 1 + h^2 + k^2; 
B = [0, (2*p/w)*sqrt(p/mu) ,0;
     sqrt(p/mu)*sin(L), sqrt(p/mu)*((w+1)*cos(L)+f)*(1/w), -sqrt(p/mu)*(h*sin(L) - k*cos(L))*(g/w);
     -sqrt(p/mu)*cos(L), sqrt(p/mu)*((w+1)*sin(L)+g)*(1/w), sqrt(p/mu)*(h*sin(L) - k*cos(L))*(f/w);
     0, 0, sqrt(p/mu)*(s_sq * cos(L))/(2*w);
     0, 0, sqrt(p/mu)*(s_sq * sin(L))/(2*w);
     0, 0, sqrt(p/mu)*(h*sin(L) - k*cos(L))*(1/w)];

A = [0 ;0 ;0 ;0 ;0 ;sqrt(mu*p)*(w/p)^2];
Xdot = [A + B*U*SI2CAN*Thr*delta/m; -Thr/c*delta];      
      
   
%% Lagrangian
Lagrangian = 0;
% Hamiltonian
H = Lagrangian + P.'*Xdot;

%% Derive costate dynamics
Pdot = -jacobian(H,X).'; % this the Euler-Lagrange equation

%%
P_1_6 = P(1:6);
Uop_num = B.'*P_1_6;
Uop_den = sqrt(Uop_num.'*Uop_num);
norm_B_lam = Uop_den;

% optimal direction
Uop = -Uop_num/Uop_den;

% construct the switching function for the Mayer-form cost functional
SF = SI2CAN*norm_B_lam*c/m + P(7);

% get the smooth engine throttle using the hyperbolic tangent smoothing (HTS) method
delta_op = 0.5*(1+tanh(SF/rho));

%% Obtain the total set of state-costate dynamics
F = [Xdot;
     Pdot];

%% Substitute for the optimal control in the RHS of the set of state-costate dynamics
F = subs(F,U,Uop);
F = subs(F,delta,delta_op);

%% State Transition Matrix
Z = [X;P];

%% This section produces a file named StCO_Dynamics_MEE.m

fc = matlabFunction(F,...
           'file','StCO_Dynamics_MEE.m',...
           'vars',{t,Z, c, Thr, rho,SI2CAN},...
           'outputs',{'Fdot'});

     

       
       
       