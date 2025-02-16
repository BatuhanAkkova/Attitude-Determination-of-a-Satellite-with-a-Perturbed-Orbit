function [r,v,nu] = Kepler_Elements_to_State(e,i,RAAN,omega,E,Period,mu,h)
% Kepler Elements to State Transformation Function

% Inputs:
% e: Eccentricity
% i: Inclination
% RAAN: Right Ascension of the Ascending Node
% omega: Argumen of Perigee
% E: Eccentric Anomaly
% Period: Orbit's Period
% mu: Earth's gravitational parameter in m^3s^-2
% h: Specific angular momentum in m^2/s

% Outputs:
% r: Position Vector in m
% v: Velocity Vector in m/s
% nu: True Anomaly

% We are using radians:
i = deg2rad(i);
RAAN = deg2rad(RAAN);
omega = deg2rad(omega);

beta = e/(1+sqrt(1-e^2)); % For true anomaly

% Preallocations below for speed:
r = zeros(3,Period);
v = zeros(3,Period);

for iter=1:Period
% Calculate true anomaly
nu(iter) = E(iter) + 2*atan(beta*sin(E(iter))/(1-beta*cos(E(iter))));

% State Vectors in perifocal frame
r_p = h^2/mu./(1+e*cos(nu(iter))).*[cos(nu(iter)) sin(nu(iter)) 0];
v_p = mu/h.*[-sin(nu(iter)) e+cos(nu(iter)) 0];

% Perifocal to ECI Transformation Matrix
R1 = [cos(-omega) -sin(-omega) 0;sin(-omega) cos(-omega) 0;0 0 1];
R2 = [1 0 0;0 cos(-i) -sin(-i);0 sin(-i) cos(-i)];
R3 = [cos(-RAAN) -sin(-RAAN) 0;sin(-RAAN) cos(-RAAN) 0;0 0 1];

% State Vector in ECI
r(:,iter) = r_p * R1 * R2 * R3;
v(:,iter) = v_p * R1 * R2 * R3;
end
r = transpose(r);
v = transpose(v);
nu = mod(nu,2*pi);