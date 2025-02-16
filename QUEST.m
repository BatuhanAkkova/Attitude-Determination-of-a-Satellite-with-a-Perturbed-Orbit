function [C_opt,q_opt] = QUEST(B_ref,B_meas,s_ref,s_meas,Period)
% QUEST Method Script
% Inputs:
% B_ref: Magnetometer output in reference frame
% B_meas: Magnetometer output in body frame
% s_ref: Sun sensor output in reference frame
% s_meas: Sun sensor output in body frame
% Period: Orbit's period

% Outputs:
% C_opt: Optimum attitude matrix 
% q_opt: Equivalent optimum quaternion

tolerance = 10e-5; % Newton-Raphson tolerance

C_opt = zeros(Period*3,3); % Preallocation for speed
q_opt = zeros(Period,4); % Preallocation for speed
P = zeros(Period*3,3); % Preallocation for speed
e = zeros(3,Period); % Preallocation for speed

ii = 1; % Matrix iteration

for i = 1:Period

% Attitude matrix calculation
vb = [transpose(B_meas(i,:)) transpose(s_meas(i,:))]; % Body frame vectors in 3xn
vi = [transpose(B_ref(i,:)) transpose(s_ref(i,:))]; % Inertial frame vectors in 3xn

weight_mag = 9.0000e-06; % var(B_meas(i,:)); % Variance of magnetometer
weight_sun = 4.0000e-06; % var(s_meas(i,:)); % Variance of sun sensor
vb = [vb(:,1)*weight_mag vb(:,2)*weight_sun]; 

B = vb*vi';

Z = [B(2,3)-B(3,2); B(3,1)-B(1,3); B(1,2)-B(2,1)];
S = B + B';
sigma = trace(B);

delta = det(S);
kappa = trace(delta*inv(S));

a = sigma^2 - kappa;
b = sigma^2 + Z'*Z;
c = delta + Z'*S*Z;
d = Z'*S^2*Z;
constant = a*b + c*sigma - d;

lambda = weight_sun + weight_mag;
last_lambda = 0.0;
while abs(lambda - last_lambda) >= tolerance
    last_lambda = lambda;
    
    f = lambda^4 - (a + b)*lambda^2 - c*lambda + constant;
    f_dot = 4*lambda^3 - 2*(a + b)*lambda - c;
    lambda = lambda - f/f_dot;
end

omega = lambda;
alpha = omega^2 - sigma^2 + kappa;
beta  = omega - sigma;
gamma = (omega + sigma)*alpha - delta;

X = (alpha*eye(3) + beta*S + S^2)*Z;
q_opt(i,:) = [X; gamma]./sqrt(gamma^2 + norm(X)^2);
q_opt(i,:) = [q_opt(i,4) q_opt(i,1) q_opt(i,2) q_opt(i,3)];
C_opt(ii:ii+2,1:3) = quat2dcm(q_opt(i,:));

% Covariance Analysis
%P(ii:ii+2,1:3) = cov(C_opt(ii:ii+2,1:3)); % The Covariance Matrix
%ee = eig(P(ii:ii+2,1:3)); % Eigenvalues of the covariance matrix
%e(1,i) = ee(1); e(2,i) = ee(2); e(3,i) = ee(3);

%ii = ii+3;
end
end