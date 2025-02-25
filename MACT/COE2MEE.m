function [p,f,g,h,k,L] = COE2MEE(e,a,i,RAAN,AoP,TA)
% Function input are the Classical Orbital Elements:
%            a - semi-major axis
%            e - eccentricity
%            i - inclination (deg)
%            RAAN - Right Ascension of the Ascending Node (deg)
%            AoP - Argument of Periapsis (deg)
%            TA - True Anomaly (deg)
% Function outputs the Modified Equinoctial Elements

% Converting degrees to radians
i = i * pi / 180;
RAAN = RAAN * pi / 180;
AoP = AoP * pi / 180;
TA = TA * pi / 180;
% Semi-Parameter - p
p = a*(1 - e^2);

% f and g
f = e*cos(AoP + RAAN);
g = e*sin(AoP + RAAN);

% h and k
h = tan(i/2)*cos(RAAN);
k = tan(i/2)*sin(RAAN);

% True Longitude (L)
L = RAAN + AoP + TA;

