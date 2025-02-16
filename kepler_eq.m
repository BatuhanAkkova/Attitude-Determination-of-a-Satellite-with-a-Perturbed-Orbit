function E = kepler_eq(e,Me,Period)
Me = mod(Me,2*pi); % Mean Anomaly within 0-2pi
E = zeros(Period,1); % Pre-allocation of Eccentric Anomaly
E(1) = Me(1) - e*sin(Me(1)) - Me(1); % First Ecc. Anomaly value
E(end) = 2*pi; % Last Ecc. Anomaly value
% Kepler Equation's Solver
for i = 1:Period-2
    E(i+1) = E(i) - (E(i) - e*sin(E(i)) - Me(i+1))/(1 - e*cos(E(i)));
end
