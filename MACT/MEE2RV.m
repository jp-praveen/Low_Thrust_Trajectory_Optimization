function [r,v] = MEE2RV(p,f,g,h,k,L,mu)
% Function input are the Modified Equinoctial elements
% Function outputs the Radius and Velocity in inertial frame

alpha_sq = h^2 - k^2;
s_sq = 1+ h^2 + k^2;
w = 1 + f*cos(L) + g*sin(L);
r_mee = p/w;

x = (r_mee/s_sq)*(cos(L) + alpha_sq*cos(L) +2*h*k*sin(L));
y = (r_mee/s_sq)*(sin(L) - alpha_sq*sin(L) +2*h*k*cos(L));
z = (2*r_mee/s_sq)*(h*sin(L) - k*cos(L));

vx = -(1/s_sq)*(sqrt(mu/p))*(sin(L) + alpha_sq*sin(L) - 2*h*k*cos(L) + g - 2*f*h*k + alpha_sq*g);
vy = -(1/s_sq)*(sqrt(mu/p))*(-cos(L) + alpha_sq*cos(L) + 2*h*k*sin(L) - f + 2*g*h*k + alpha_sq*f);
vz = (2/s_sq)*(sqrt(mu/p))*(h*cos(L) + k*sin(L) + f*h + g*k);

r = [x;y;z];
v = [vx;vy;vz];