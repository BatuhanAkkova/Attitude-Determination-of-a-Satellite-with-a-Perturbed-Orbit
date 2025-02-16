function a_srp = solar_p(m,A,r,jd)
% Solar-Radiation Pressure Perturbation
% Inputs:
% m: Mass of the Satellite
% A: Cross-Sectional Area of the Satellite
% r: Position Vector
% jd: Julian Dates

% Output: Acceleration due to the Solar Radiation Pressure in m/s^2

c_r = 1; % Reflectivity Constant, assume 1
p_srp = 4.57E-6; % Solar Pressure per unit area in N/m2

r_sun = planetEphemeris(jd,'Earth','Sun'); % Sun Position from Earth

r_sc_sun = r_sun - r; % Sun Position from Spacecraft

a_srp = -p_srp*c_r*A/m*r_sc_sun./vecnorm(r_sc_sun,2,2);