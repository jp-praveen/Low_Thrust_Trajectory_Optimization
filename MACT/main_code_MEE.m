% Author: Praveen
% This code aims at solving heliocentric fuel-optimal trajectories 
% from the Earth to Dionysus using three-dimensional dynamics

clear 
clc
close all

tic
%%
%==========================================================================
%Part - A: Mapping of Cartesian Costates to MEE costate using Symbolic vars
%==========================================================================
mu = 132712440018; % Sun gravitational constant (km^3/s^2)
AU = 149.6e6;      % One astronomical unit (km)
DU = AU;                      % Distance unit
TU = DU^1.5/sqrt(mu);         % Time unit
VU = DU/TU;                   % Velocity unit

syms x y z vx vy vz r v rmag vmag E a e emag RAAN AoP TA h hmag N Nmag inc p f g hmee L v_radial Q parQ 
mu = 1;

% Defining Orbital Elements in terms of Cartesian Coordinates
r = [x; y; z];
v = [vx; vy; vz];
rmag = norm(r);
vmag = norm(v);

% Total Energy 'E'
E = vmag^2/2 - mu/rmag;

% Semi-Major Axis ()
a = -mu/(2*E);

% Eccentricity Vector
e = (1/mu)*((vmag^2 - (mu/rmag))*r - dot(r,v)*v);
emag = norm(e);

% Specific Angular Momentum
h = cross(r,v);
%hmag = norm(h);
hmag = sqrt(h(1)^2 + h(2)^2 + h(3)^2);  

% Radial velocity
v_radial = dot(r,v)/rmag; 

% Node Vector
khat = [0;0;1];
N = cross(khat,h);
Nmag = norm(N);

% RAAN in radians
%RAAN = acos(N(1)/Nmag);
Nsubs = subs(N,{x y z vx vy vz},{-3637871.08165087/DU 147099798.784179/DU -2261.44104198769/DU -30.2650979882182/VU -0.848685467901138/VU 5.05303606281563e-05/VU});

% Choosing the correct quadrant for RAAN in degrees
if Nsubs(2) < 0
    RAAN = 2*pi - acos(N(1)/Nmag);
else
    RAAN = acos(N(1)/Nmag);
end

% Argument of Periapis in radians
%AoP = acos(dot(N,e)/(Nmag*emag));
esubs = subs(e,{x y z vx vy vz},{-3637871.08165087/DU 147099798.784179/DU -2261.44104198769/DU -30.2650979882182/VU -0.848685467901138/VU 5.05303606281563e-05/VU});

% Choosing the correct quadrant for Argument of Periapsis
if esubs(3) < 0
    AoP = 2*pi - acos(dot(N,e)/(Nmag*emag));
else 
    AoP = acos(dot(N,e)/(Nmag*emag));
end

% Inclination in radians
inc = acos(h(3)/hmag);

% True Anomaly in radians
%TA = acos(dot(e,r)/(emag*rmag));
v_radialsubs = subs(v_radial,{x y z vx vy vz},{-3637871.08165087/DU 147099798.784179/DU -2261.44104198769/DU -30.2650979882182/VU -0.848685467901138/VU 5.05303606281563e-05/VU});

if v_radialsubs < 0
    TA = 2*pi - acos(dot(e,r)/(emag*rmag));
else
    TA = acos(dot(e,r)/(emag*rmag));
end

% Defining MEE in terms of Orbital Elements ()
p = a*(1 - emag^2);

f = emag*cos(RAAN + AoP);

g = emag*sin(RAAN + AoP);

hmee = tan(inc/2)*cos(RAAN);

k = tan(inc/2)*sin(RAAN); 

L = RAAN + AoP + TA;

% Q matrix
Q = [p; f; g; hmee; k; L];

% Q Derivative
parQ = [diff(Q,x),diff(Q,y),diff(Q,z),diff(Q,vx),diff(Q,vy),diff(Q,vz)];

parQsubs = subs(parQ,{x y z vx vy vz},{-3637871.08165087/DU 147099798.784179/DU -2261.44104198769/DU -30.2650979882182/VU -0.848685467901138/VU 5.05303606281563e-05/VU});
parQeval = eval(parQsubs);

%%
%==========================================================================
%Part - B: Mapping guessed Cartesian costates using ACT method to MEE
%costates 
%==========================================================================

clearvars -except parQeval DU TU VU
% Guess the Cart costates

%angle_i = 1*rand(7,1);
%angle_i = [1.542103265572267e3;-2.657072594277044e2;0;0;0;0;0.0];
angle_i = [1000*rand;10*rand;0;0;0;0;0.5];

lvmag = angle_i(1);
lvmag_dot  = angle_i(2); 
alp  = angle_i(3);
alp_dot =  angle_i(4);
beta = angle_i(5);
beta_dot = angle_i(6);
lm = angle_i(7);


x = -3637871.08165087/DU;             % km
y = 147099798.784179/DU;              % km
z = -2261.44104198769/DU;             % km

vx = -30.2650979882182/VU;              % km/s
vy = -0.848685467901138/VU;             % km/s
vz = 5.05303606281563e-05/VU;           % km/s

% Velocity vector
Vvec = [vx;vy;vz];

% Radius vector
Rvec = [x;y;z];

% Angular momentum vector
hvec = cross(Rvec,Vvec);

% 
Vdot = [-x/(x^2 + y^2 + z^2)^1.5; 
        -y/(x^2 + y^2 + z^2)^1.5; 
        -z/(x^2 + y^2 + z^2)^1.5]; 
   
Rdot = [vx;
        vy;
        vz];

hdot = cross(Rdot,Vvec) + cross(Rvec,Vdot);

% Spacecraft centered frame's unit vectors
Vhat = (Vvec/norm(Vvec));
hhat = (hvec/norm(hvec));
bhat = cross(hhat,Vhat);

Vhat_dot = Vdot/norm(Vvec) - Vvec*(dot(Vvec,Vdot)/norm(Vvec))/norm(Vvec)^2;

hhat_dot = hdot/norm(hvec) - hvec*(dot(hvec,hdot)/norm(hvec))/norm(hvec)^2;

bhat_dot = cross(hhat_dot,Vhat) + cross(hhat,Vhat_dot);


% Thrust pointing vector in Spacecraft centered frame
Utv = [cos(alp)*cos(beta);sin(alp)*cos(beta);sin(beta)];
Utv_dot = [-alp_dot*sin(alp)*cos(beta)-beta_dot*cos(alp)*sin(beta); 
            alp_dot*cos(alp)*cos(beta)-beta_dot*sin(alp)*sin(beta); 
            beta_dot*cos(beta)];

% Direction Cosine Matrix (DCM) to transform from Spacecraft centered
% frame to Cartesian frame ()
D = [dot([1;0;0],Vhat), dot([1;0;0],bhat), dot([1;0;0],hhat);
     dot([0;1;0],Vhat), dot([0;1;0],bhat), dot([0;1;0],hhat);
     dot([0;0;1],Vhat), dot([0;0;1],bhat), dot([0;0;1],hhat)]; 
 
Ddot = [dot([1;0;0],Vhat_dot), dot([1;0;0],bhat_dot), dot([1;0;0],hhat_dot);
        dot([0;1;0],Vhat_dot), dot([0;1;0],bhat_dot), dot([0;1;0],hhat_dot);
        dot([0;0;1],Vhat_dot), dot([0;0;1],bhat_dot), dot([0;0;1],hhat_dot)];

% Thrust pointing vector in Cartesian frame    
Utc = D*Utv;
Utc_dot = Ddot*Utv + D*Utv_dot;

% Velocity and Position costates
lv = -lvmag*Utc;
lvmag_dot = -1/(Utc.'*Vvec) * (lvmag*Utc_dot.'*Vvec - lv.'*Vdot);
lvdot = -lvmag_dot*Utc - lvmag*Utc_dot;

lr = -lvdot;

lambda_cart = [lr;lv;lm];

parQeval_trans = transpose(parQeval);
lambda_MEE = inv(parQeval_trans)*lambda_cart(1:6);
lambda_MEE = [lambda_MEE;lm];

%%
%==========================================================================
%Part - C: Solving TPBVP using the mapped MEE costates 
%==========================================================================
clearvars -except lambda_MEE lambda_cart angle_i  
% define constants

mu = 132712440018; % Sun gravitational constant (km^3/s^2)
AU = 149.6e6;      % One astronomical unit (km)
Isp = 3000;       % seconds
g0 = 9.80665;     % m/s^2
T = 0.32;         % Max thrust in N

% defining Canonical units for scaling 

DU = AU;                      % Distance unit

TU = DU^1.5/sqrt(mu);         % Time unit

VU = DU/TU;                   % Velocity unit

eta = TU^2/(DU*1000);      % conversion factor from m/s^2 to DU/TU^2

c = Isp * g0 / TU;         % scaling c

mu = 1;

% Initial states

R_i =[-3637871.08165087;
      147099798.784179;
      -2261.44104198769] ;  % km

x_i = R_i(1);      % km
y_i = R_i(2);      % km
z_i = R_i(3);      % km

v_xi = -30.2650979882182;      % km/s
v_yi = -0.848685467901138;     % km/s
v_zi = 5.05303606281563e-05;   % km/s

% Cartesian to Classical Orbital elements for the initial states
pos_i = [x_i;y_i;z_i];
vel_i = [v_xi;v_yi;v_zi];
[e_i,a_i,i_i,RAAN_i,AoP_i,TA_i] = RV2COE(pos_i/DU,vel_i/VU,1);

% Classical elements to Equinoctial elements for the initial states
[p_i,f_i,g_i,h_i,k_i,L_i] = COE2MEE(e_i,a_i,i_i,RAAN_i,AoP_i,TA_i);

m_i = 4000;                   % kg
t_i = 0;

% Scaled initial state vector
S_i = [p_i;f_i;g_i;h_i;k_i;L_i;m_i];   

% Final states

x_f = -302452014.884247;             % km
y_f = 316097179.632028;              % km
z_f = 82872290.075518;                % km

v_xf = -4.53347379984029;            % km/s
v_yf = -13.1103098008475;            % km/s
v_zf = 0.65616382601745;            % km/s

% Cartesian to Classical Orbital elements for the final states
pos_f = [x_f;y_f;z_f];
vel_f = [v_xf;v_yf;v_zf];
[e_f,a_f,i_f,RAAN_f,AoP_f,TA_f] = RV2COE(pos_f/DU,vel_f/VU,1);

% Classical elements to Equinoctial elements for the final states
[p_f,f_f,g_f,h_f,k_f,L_f] = COE2MEE(e_f,a_f,i_f,RAAN_f,AoP_f,TA_f);

t_f = 3534 * 86400;        % time of flight in seconds

t_f = t_f/TU;                 % scaled time of flight

% Final Longitude must be greater than the initial Longitude
while L_f<L_i
    L_f = L_f + 2*pi;    
end

% Defining number of revolutions
N_rev = 5;
L_f = L_f + 2*pi*N_rev;

% Scaled final state vector
S_f = [p_f;f_f;g_f;h_f;k_f;L_f];        

% propagation time interval
tspan=[0,t_f];

% set solver options
options_ode    = odeset('RelTol',1e-10,'AbsTol',1e-10);
options_fsolve = optimoptions('fsolve','Display','iter','SpecifyObjectiveGradient',false,'MaxFunEvals',500000,'MaxIter',500000,'TolFun',1e-14,'TolX',1e-14,'Algorithm','Levenberg-Marquardt'); % Option to display output

% Data Stucture to be passed to functions
rho = 1;
PD.rho = rho;
PD.c  =  c;
PD.T = T;     
PD.tspan = tspan;
PD.eta = eta;
PD.S_i = S_i;
PD.S_f = S_f;
PD.options_ode = options_ode;

% Initiating the norm of residual
norm_res=1;

flag = 1;
mass_arr = [];
rho_arr = [];
while norm_res>1e-5
    flag = flag + 1;
    CS_i = lambda_MEE;
    [x_c,fval,exitflag,output] = fsolve(@BC_jaco_MEE,CS_i,options_fsolve,PD);
    %[x_c,fval] = fsolve(@BC_jaco_MEE,CS_i,PD);
    norm_res = sqrt(fval.'*fval);
    if flag > 1
        break
    end
end

%%
rho_arr = [rho_arr rho];
CS_i = x_c;
X_i = [ S_i;   % states
       CS_i;];  % costates

[time,X_p] = ode45(@StCO_Dynamics_MEE,tspan,X_i,options_ode,c, T, rho,eta); 

% Storing the final mass value for each rho value for plotting
mass_arr = [mass_arr X_p(end,7)];

norm_res = 1;
rho = rho*0.1;
k1 = 1;

% Reducing the rho value iteratively to find the solution of the 
% original problem. 
while norm_res>1e-5
    
    PD.rho = rho;
    [x_c,fval,exitflag,output] = fsolve(@BC_jaco_MEE,CS_i,options_fsolve,PD);
    norm_res = sqrt(fval.'*fval);
    
    % Checking for convergance
    if norm_res > 1e-5
       
        rho = rho*1.001;
        continue
    end
    
    % Converged solution for a rho value being used as the initial guess
    % for the next reduced rho value
    CS_i = x_c;
    
    % Reducing the value of rho 
    rho = rho*0.1;
    
    % Reinitializing the norm of residual
    norm_res = 1;
    
    if rho < 1e-3
        break
    end
    
    %Storing the rho value for plotting purposes
    rho_arr = [rho_arr rho];
    
    CS_i = x_c;
    X_i = [ S_i;   % states
            CS_i;];  % costates

    [time,X_p] = ode45(@StCO_Dynamics_MEE,tspan,X_i,options_ode,c, T, rho,eta); 
    
    % Storing the final mass value for each rho value for plotting
    mass_arr = [mass_arr X_p(end,7)];
end
%%
rho = rho/0.1;
X_i = [ S_i;   % states
       CS_i;];  % costates

[time,X_p] = ode45(@StCO_Dynamics_MEE,tspan,X_i,options_ode,c, T, rho,eta);    
   
%% Extracting the final solution

p_f  = X_p(:,1);
f_f  = X_p(:,2);
g_f  = X_p(:,3);
h_f = X_p(:,4);
k_f = X_p(:,5);
L_f = X_p(:,6);
m_f = X_p(:,7);
lp_f  = X_p(:,8);
lf_f  = X_p(:,9);
lg_f  = X_p(:,10);
lh_f = X_p(:,11);
lk_f = X_p(:,12);
lL_f = X_p(:,13);
lm_f = X_p(:,14);

x = [];
y = [];
z = [];
vx = [];
vy = [];
vz = [];
r = [];

% Converting MEE to Radius and Velocity vectore
for k1 = 1:length(p_f)
    [pos,vel] = MEE2RV(p_f(k1),f_f(k1),g_f(k1),h_f(k1),k_f(k1),L_f(k1),mu);
    x(k1) = pos(1)*DU;
    y(k1) = pos(2)*DU;
    z(k1) = pos(3)*DU;
    vx(k1) = vel(1)*VU;
    vy(k1) = vel(2)*VU;
    vz(k1) = vel(3)*VU;
    r_vec(:,k1) = [x(k1);y(k1);z(k1)];
    v_vec(:,k1) = [vx(k1);vy(k1);vz(k1)];
    r(k1) = sqrt(x(k1)^2 + y(k1)^2 + z(k1)^2);
    v(k1) = sqrt(vx(k1)^2 + vy(k1)^2 + vz(k1)^2);
    
    % Direction Unit Vectors
    ir(:,k1) = r_vec(:,k1)/r(k1);
    in(:,k1) = cross(r_vec(:,k1),v_vec(:,k1))/norm(cross(r_vec(:,k1),v_vec(:,k1)));
    it(:,k1) = cross(in(:,k1),ir(:,k1));
end

% Finding the optimal  control
for k1 = 1:length(p_f)
    w_f(k1) = 1 + f_f(k1)*cos(L_f(k1)) + g_f(k1)*sin(L_f(k1));  
    s_sq_f(k1) = 1 + h_f(k1)^2 + k_f(k1)^2; 
    B(:,:,k1) = [0, (2*p_f(k1)/w_f(k1))*sqrt(p_f(k1)/mu) ,0;
     sqrt(p_f(k1)/mu)*sin(L_f(k1)), sqrt(p_f(k1)/mu)*((w_f(k1)+1)*cos(L_f(k1))+f_f(k1))*(1/w_f(k1)), -sqrt(p_f(k1)/mu)*(h_f(k1)*sin(L_f(k1)) - k_f(k1)*cos(L_f(k1)))*(g_f(k1)/w_f(k1));
     -sqrt(p_f(k1)/mu)*cos(L_f(k1)), sqrt(p_f(k1)/mu)*((w_f(k1)+1)*sin(L_f(k1))+g_f(k1))*(1/w_f(k1)), sqrt(p_f(k1)/mu)*(h_f(k1)*sin(L_f(k1)) - k_f(k1)*cos(L_f(k1)))*(f_f(k1)/w_f(k1));
     0, 0, sqrt(p_f(k1)/mu)*(s_sq_f(k1) * cos(L_f(k1)))/(2*w_f(k1));
     0, 0, sqrt(p_f(k1)/mu)*(s_sq_f(k1) * sin(L_f(k1)))/(2*w_f(k1));
     0, 0, sqrt(p_f(k1)/mu)*(h_f(k1)*sin(L_f(k1)) - k_f(k1)*cos(L_f(k1)))*(1/w_f(k1))];
end

for k2 = 1:length(p_f)
    B_temp = B(:,:,k2);
    Uop_num(:,k2) = B_temp.'*[lp_f(k2);lf_f(k2);lg_f(k2);lh_f(k2);lk_f(k2);lL_f(k2)];
    norm_B_lam(k2) = sqrt(Uop_num(:,k2).'*Uop_num(:,k2));
    Uop_LF(:,k2) = -Uop_num(:,k2)/norm_B_lam(k2);
    Q = [ir(:,k2),it(:,k2),in(:,k2)];
    u_cart(:,k2) = Uop_LF(:,k2);
    Uop(:,k2) = Q*Uop_LF(:,k2);
    ur(k2) = Uop(1,k2);
    ut(k2) = Uop(2,k2);
    un(k2) = Uop(3,k2);
    SF(k2) = eta*norm_B_lam(k2)*c/m_f(k2) + lm_f(k2);
    delta(k2) = 0.5*(1+tanh(SF(k2)/rho));
    alpha(k2) = atan2(ut(k2),ur(k2));
    beta(k2) = atan2(un(k2),sqrt(ur(k2)^2 + ut(k2)^2));
end

time = time * TU / 86400;

%% Plotting

% r and v plot
figure(1);
subplot(2,1,1);
plot(time,r,'r');
xlabel('Time(days)')
ylabel('r (km)')
subplot(2,1,2);
plot(time,v,'g');
xlabel('Time(days)')
ylabel('v (km/s)')

% Cartesian position coordinates plot
figure(2);
subplot(3,1,1);
plot(time,x,'r');
xlabel('Time(days)')
ylabel('x (km)')
subplot(3,1,2);
plot(time,y,'g');
xlabel('Time(days)')
ylabel('y (km) ')
subplot(3,1,3);
plot(time,z,'b');
xlabel('Time(days)')
ylabel('z (km) ')

% Plot of Final mass
figure(3);
plot(time,m_f,'b');
xlabel('Time(days)')
ylabel('mass (kg)')

% Switching function and delta plot
figure(4);
plot(time,T*delta,'r','LineWidth',2);
yyaxis left
xlabel('Time(days)')
ylabel('Thrust')
yyaxis right
plot(time,SF,'g','LineWidth',2);
%xlabel('Time(days)')
ylabel('Switching Function')
legend('Thrust','Switching Function')

% Position of Earth and Dionysus
x_earth = x(1);
y_earth = y(1);
z_earth = z(1);
x_dionysus = x(end);
y_dionysus = y(end);
z_dionysus = z(end);

% Velocity of Earth and Dionysus
vx_earth = vx(1);
vy_earth = vy(1);
vz_earth = vz(1);
vx_dionysus = vx(end);
vy_dionysus = vy(end);
vz_dionysus = vz(end);

% State vector of Earth and Dionysus

X_earth = [x_earth;y_earth;z_earth;vx_earth;vy_earth;vz_earth];
X_dionysus = [x_dionysus;y_dionysus;z_dionysus;vx_dionysus;vy_dionysus;vz_dionysus];
tspan_earth = [0,365*86400];
tspan_dionysus = [0,1200*86400];

options_ode    = odeset('RelTol',1e-10,'AbsTol',1e-10);
[time_earth,X_p_earth] = ode45(@planet_orbit,tspan_earth,X_earth,options_ode); 
[time_dionysus,X_p_mars] = ode45(@planet_orbit,tspan_dionysus,X_dionysus,options_ode);

% Trajectory Plot
figure(5);
plot3(x/AU,y/AU,z/AU,'k','LineWidth',2);
xlabel('x(AU)')
ylabel('y(AU)')
zlabel('z(AU)')
hold on
plot3(X_p_earth(:,1)/AU,X_p_earth(:,2)/AU,X_p_earth(:,3)/AU,'-c','LineWidth',2)
plot3(X_p_mars(:,1)/AU,X_p_mars(:,2)/AU,X_p_mars(:,3)/AU,'-g','LineWidth',2)
n=5;
x_quiver = x(1:n:end)/AU;
y_quiver = y(1:n:end)/AU;
z_quiver = z(1:n:end)/AU;

% Control Profile
ux = u_cart(1,:);
uy = u_cart(2,:);
uz = u_cart(3,:);

% Unit thrust steering vector
ux_quiver = ur(1:n:end)/AU;
uy_quiver = ut(1:n:end)/AU;
uz_quiver = un(1:n:end)/AU;
delta_quiver = delta(1:n:end);
quiver3(x_quiver,y_quiver,z_quiver,ux_quiver.*delta_quiver,uy_quiver.*delta_quiver,uz_quiver.*delta_quiver,0.5,'b','LineWidth',1)
plot3(x(1)/AU,y(1)/AU,z(1)/AU,'.','MarkerEdgeColor','r','MarkerFaceColor',[.49 1 .63],'MarkerSize',20)
plot3(x(end)/AU,y(end)/AU,z(end)/AU,'.','MarkerEdgeColor','m','MarkerFaceColor',[.49 1 .63],'MarkerSize',20)
legend('Trajectry','Earth Orbit', 'Dionysus Orbit','Thrust direction','Departure Point','Arrival Point','Location','northeast');

% Final Mass Vs rho
figure(6)
loglog(rho_arr)
plot(rho_arr,mass_arr,'.-')
set(gca, 'XScale', 'log')
xlabel('\rho')
ylabel('m_f(kg)')
hold on

end_time = toc



