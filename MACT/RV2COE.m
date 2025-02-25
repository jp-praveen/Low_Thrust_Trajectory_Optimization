function [e_mag,a,i,RAAN,AoP,TA] = RV2COE(pos,vel,mu)
% Function input are the position and velocity column vectors, and
% gravitational parameter
% Function outputs the classical orbital elements
 

% Position in cartesian  form
x = pos(1);
y = pos(2);
z = pos(3);
r = pos;

% Velocity in cartesian form
vx = vel(1);
vy = vel(2);
vz = vel(3);
v = vel;

% 
r_mag = sqrt(x^2 + y^2 + z^2);
v_mag = sqrt(vx^2 + vy^2 + vz^2);

% Radial velocity
v_radial = dot(r,v)/r_mag; 

% Specific Angular Momentum
h = cross(pos,vel);
h_mag = sqrt(h(1)^2 + h(2)^2 + h(3)^2);  

% Inclination in degrees
i = acosd(h(3)/h_mag); % an orbital element! (1)

% Node Vector
k_hat = [0;0;1];
N = cross(k_hat,h);
N_mag = sqrt(N(1)^2 + N(2)^2 + N(3)^2);

% RAAN in degrees
RAAN = acosd(N(1)/N_mag); % an orbital element! (2)

% Choosing the correct quadrant for RAAN in degrees
if N(2) < 0
    RAAN = 360 - RAAN;
else
    RAAN = RAAN;
end

% Eccentricity vector
e = (1/mu)*((v_mag^2 - (mu/r_mag))*r - dot(r,v)*v);
e_mag = sqrt(e(1)^2 + e(2)^2 + e(3)^2); % an orbital element! (3)

% Argument of Periapis in degrees
AoP = acosd(dot(N,e)/(N_mag*e_mag)); % an orbital element! (4)

% Choosing the correct quadrant for Argument of Periapsis
if e(3) < 0
    AoP = 360 - AoP;
else 
    AoP = AoP;
end

% Semi-Major axis in Km
E_tot = v_mag^2/2 - mu/r_mag; % Total Orbital energy

if e_mag ~= 1 
    a = -mu/(2*E_tot);        % an orbital element! (5)
else
    a = inf;
end

% True Anomaly in degrees
TA = acosd(dot(e,r)/(e_mag*r_mag));   % an orbital element! (6)

if v_radial < 0
    TA = 360 - TA;
end
    


    

   
















