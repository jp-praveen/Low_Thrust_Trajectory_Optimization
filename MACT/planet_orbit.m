function dstate_vector = planet_orbit(t,X)

mu = 1.327e11;  % Sun Gravitational parameter in kg^3/m^2

% Initial radius and velocity

x  = X(1);
y  = X(2); 
z  = X(3);
vx =  X(4);
vy = X(5);
vz = X(6);

r = sqrt(x^2+y^2+z^2);

% Position Differential equations
dx = vx;
dy = vy;
dz = vz;

% Velocity Differential equations
dvx = -mu*x/r^3;
dvy = -mu*y/r^3;
dvz = -mu*z/r^3;

dstate_vector = [dx; dy; dz; dvx; dvy; dvz]; 