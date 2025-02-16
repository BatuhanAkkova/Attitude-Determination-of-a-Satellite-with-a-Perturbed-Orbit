function y = Encke(r,v,mu,h,Period,a_p)
% Encke's Method
% Inputs: 
% r: Position Vector
% v: Velocity Vector
% mu: Gravitational parameter of Earth
% h: Specific Angular Momentum
% Period: Orbit's Period
% a_p: Acceleration due to Perturbations

% Output: State Vectors (both r and v)

t0 = 0; tf = Period; % Simulation initial and final time in s

% Initial State Vector
r0 = r(1,:);
v0 = v(1,:);

del_t = 1; % Step Size
opt = odeset('MaxStep',del_t);

t = t0; % Initial Time for iter.
tsave = t0; % Initial Soln. Time Vector
y = [r0 v0]; % Initial state vector for iter.
del_y0 = zeros(6,1); % Initial perturbation vector for iter.

t = t + del_t; % First Time Step
i = 1;

while t <=tf + del_t/2 

% Integrate rates from t0 to t
[~,z] = ode45(@rates,[t0 t],del_y0,opt);

% Compute osculating state
[rosc,vosc] = State(r0,v0,mu,t-t0,h);

% Rectifying
r0 = rosc + z(end,1:3);
v0 = vosc + z(end,4:6);
t0 = t;

% Next Step
tsave = [tsave;t];
y = [y; [r0 v0]];
t = t + del_t;
del_y0 = zeros(6,1);
i = i + 1;
end
t = tsave;

function dfdt = rates(t,f)
del_r = f(1:3)'; % Position Deviation
del_v = f(4:6)'; % Velocity Deviation

% Osculating State
[rosc,vosc] = State(r0,v0,mu,t-t0,h);

% Perturbed State
rpp = rosc + del_r;
vpp = vosc + del_v;

% Magnitudes
Rosc = norm(rosc);
Rpp = norm(rpp);

% Total Acc.
F = 1 - (Rosc/Rpp)^3;
del_a = -mu/Rosc^3*(del_r - F*rpp) + a_p(i,:);
dfdt = [del_v(1) del_v(2) del_v(3) del_a(1) del_a(2) del_a(3)]';
end
end
