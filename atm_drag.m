% Atmospheric Drag Perturbation
function adrag = atm_drag(rho,v,m,A)
% Inputs:
% rho = density at given altitude
% v = state velocity of satellite
% m: Mass of the Satellite
% A: Cross-Sectional Area of the Satelite

% Output: Acceleration due to the Atmospheric Drag in m/s^2
Cd = 0.88; % Drag Coefficient Estimated from Research

v(:,1:2) = v(:,1:2) - 465.1; % satellite vel. + Earth rot. velocity in m/s
v_mag = vecnorm(v,2,2); % magnitude of v vector

adrag = -1/2.*rho*(Cd*A/m).*v_mag.*v;
end




