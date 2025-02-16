function q_dot = EOM(I,L,w0,q0,Period)
% Equations of Motion Script
% Inputs: 
% I: Moment of inertia
% L: Control Torque
% w0: Initial angular velocity
% q0: Initial quaternion
% Period: Orbit's period

% Output: propagated quaternion

tspan = 1:Period; % Step Size

% Dynamic Equation of Motion
% derivative of w from the book Analytical Mechanics of Space Systems by
% Junkins and Schaub
w_dot_f = @(t,w) [(-(I(3)-I(2))*w(2)*w(3)+L)/I(1);... 
    (-(I(1)-I(3))*w(3)*w(1)+L)/I(2);... 
    (-(I(2)-I(1))*w(1)*w(2)+L)/I(3)];

% Solution to w dot functions with 4th Order Runge-Kutta Method
options = odeset('MaxStep',1);
w_dot = ode45(w_dot_f, tspan, w0, options);
w_dot = transpose(w_dot.y);

% Kinematic Equation of Motion
% The q dot function from book Analytical Mechanics of Space Systems by
% Junkins and Schaub
q_dot_f =@(t,q) 1/2.*[0 -w_dot(1) -w_dot(2) -w_dot(3); w_dot(1) 0 w_dot(3) -w_dot(2); w_dot(2) -w_dot(3) 0 w_dot(1); w_dot(3) w_dot(2) -w_dot(1) 0]*q;

% Solving the q dot function with 4th Order Runge-Kutta Method:
options = odeset('MaxStep',1);
q_dot = ode45(q_dot_f, tspan, q0, options);
q_dot = transpose(q_dot.y);
end
