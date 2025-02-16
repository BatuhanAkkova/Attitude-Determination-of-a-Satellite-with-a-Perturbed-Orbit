function [e_pert,inc_pert,RAAN_pert,omega_pert,h_pert,nu_pert] = State_to_Kepler_Elements(r,v,mu,Period)
% State Vector to Kepler Elements Transformation Function
% Inputs:
% r: Position Vector
% v: Velocity Vector
% mu: Gravitational parameter of Earth
% Period: Orbit's Period

% Outputs: The Kepler Parameters

% Preallocations
h_pert = zeros(Period,1);
inc_pert = zeros(Period,1);
RAAN_pert = zeros(Period,1);
e_pert = zeros(Period,1);
omega_pert = zeros(Period,1);
nu_pert = zeros(Period,1);

for iter = 1:Period

r_mag = norm(r(iter,:)); % Magnitude of Position Vector

v_r = dot(r(iter,:)/r_mag, v(iter,:)); % Radial Velocity

h_vec = cross(r(iter,:), v(iter,:)); % Angular Momentum Vector
h_pert(iter) = norm(h_vec); % Angular Momentum

inc_pert(iter) = acos(h_vec(3) / h_pert(iter)); % Inclination

K = [0 0 1]; % Z-Axis
N_vec = cross(K,h_vec); % Ascending Node Vector
N = norm(N_vec); % Ascending Node

if N_vec(2) < 0 % If y - component is neg.
RAAN_pert(iter) = 2*pi - acos(N_vec(1)/N); % Right asc..
else
RAAN_pert(iter) = acos(N_vec(1)/N);
end

e_vec = cross(v(iter,:), h_vec) / mu - r(iter,:)/r_mag; % Eccentricity Vector
e_pert(iter) = norm(e_vec); % Eccentricity

if e_vec(3) < 0 % If z - component is neg.
omega_pert(iter) = 2*pi - acos(dot(N_vec, e_vec) / (N * e_pert(iter))); % Arg. of Perigee
else
omega_pert(iter) = acos(dot(N_vec, e_vec) / (N * e_pert(iter)));
end

nu_pert = 0:2*pi/Period:2*pi;
nu_pert = nu_pert(1:Period); % True Anomaly

end
