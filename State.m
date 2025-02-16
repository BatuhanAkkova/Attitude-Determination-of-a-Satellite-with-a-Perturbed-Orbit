function [r,v] = State(r0,v0,mu,dtrue,h) 
% State Vector Propagation From Initial State
% Inputs:
% r0: Position Vector
% v0: Velocity Vector
% mu: Gravitational Parameter of Earth
% dtrue: Time difference
% h: Specific Angular Momentum

% Outputs:
% r: Position Vector
% v: Velocity Vector

[f,g,fdot,gdot] = Lagrange_Coeff(mu,r0,v0,h,dtrue);

r = f*r0 + g*v0;
v = fdot*r0 + gdot*v0;

end