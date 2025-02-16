function [f,g,fdot,gdot] = Lagrange_Coeff(mu,r0,v0,h,dtrue)
% Lagrange Coefficients and their Derivatives Calculation
% Inputs:
% mu: Gravitational Parameter of Earth
% r0: Position Vector
% v0: Velocity Vector
% h: Specific Angular Momentum
% dtrue: Time difference

% Outputs:
% f and g: Lagrange Coefficients
% fdot and gdot: Their Derivatives

vr0 = dot(v0,r0)/norm(r0); % Radial v0
r0 = norm(r0); % Mag. of r0
s = sind(dtrue); % Sine of dt
c = cosd(dtrue); % Cosine of dt

r = h^2/mu/(1 + (h^2/mu/r0 - 1)*c - h*vr0*s/mu);

f = 1 - mu*r*(1 - c)/h^2;
g = r*r0*s/h;

fdot = mu/h*(vr0/h*(1 - c) - s/r0);
gdot = 1 - mu*r0/h^2*(1 - c);

end