% MY INITIAL VALUES

% Earth's Values
Re = 6356752.3; % Earth's radius in m
rot_E = 7.292115*1e-5; % Earth's rotational angular velocity as rad/s
mu = 3.986004418*10^14; % Earth's gravitational parameter in m^3s^-2

% Orbital Values
semi_major = ((Re + 300E3) + (Re + 2000E3))/2; % semi-major axis in m, 300km rp 2000km ra
e = 1 - (Re+300E3)/semi_major; % eccentricity
inc = 45; % inclination in deg
RAAN = 60; % right ascension of the ascending node in deg
omega = 30; % argument of perigee in deg
Period = fix(2*pi*semi_major^(3/2)/sqrt(mu)); % in seconds
dt = transpose(1:Period);
Me = 2*pi/Period*dt; % in radians
h = sqrt(mu*semi_major*(1-e^2)); % Specific angular momentum in m^2/s
dtrue = 0.1; % Change in True Anomaly for Iterations in deg

% Satellite Initial Values
L = 3.6E-10; % Torque Vector
I = [0.0021;0.002;0.0019]; % Inertia Matrix (Diag)
w0 = [0.0011 0.0012 0.0013]; % Initial angular velocities for ode
q0 = [0 0.018 0.009 0.045]; % Initial quaternion for ode
q0(1) = sqrt(1 - q0(2)^2 - q0(3)^2 - q0(4)^2); % Scalar part of the quaternion 

% Selected magnetometer: https://www.satnow.com/products/magnetometers/antrix-corporation-limited/40-1183-digital-miniature-magnetometer
noise_mag = randn; % Normally distributed scalar number for creating noise
deviation_mag = 0.003; % Non-linearity of our magnetometer

% Selected sun sensor: https://www.aac-clyde.space/what-we-do/space-products-components/adcs/ss200
noise_sun = randn; % Normally distributed scalar number for creating noise
deviation_sun = 0.002;

weight = 1/2; % For FOAM Method Script

% Non-Spherical Earth Coefficients
J2 = 1082.6267/10^6;
J3 = -2.5327/10^6;
J4 = -1.6196/10^6;
J5 = -0.2273/10^6;
J6 = -0.6083/10^8;

% Satellite Physical, 1kg 10x10x10 cm
m = 1;
A = 0.1*0.1;
