% Third-Body Perturbations
function [a_moon] = Third_Body(r,jd,Period)
% Inputs:
% r: State Position
% jd: Julian Dates
% Period: Orbit's Period

% Output: Acceleration due to the Moon in m/s^2

mu_moon = 4.9048695*10^12; % Gravitational parameter of Moon in m^3/s^2

% Preallocation for speed
a_moon = zeros(Period,3);
alpha = zeros(Period,1);
r_s3 = zeros(Period,3);
r_s3_mag = zeros(Period,1);

r_3 = (planetEphemeris(jd,'Earth','Moon'))*10^3; % position of Moon to Earth in m
r_3_mag = vecnorm(r_3,2,2); % magnitude of the above

r_mag = vecnorm(r,2,2); % magnitude of the position of satellite to Earth

for i=1:Period

alpha(i,1) = acos(dot(r(i,:),r_3(i,:))./(norm(r(i,:)).*norm(r_3(i,:)))); % alpha = acos((r*r3)/(r_mag*r3_mag))

r_s3(i,:) = sqrt(abs(r(i,:).^2 + r_3(i,:).^2 - 2.*r_mag(i,:).*r_3_mag(i,:).*cos(alpha(i,1)))); % The position of Moon with respect to our spacecraft
r_s3_mag(i,:) = vecnorm(r_s3(i,:),2,2);

a_moon(i,:) = mu_moon*(r_s3(i,:)./(r_s3_mag(i,:)).^3 - r_3(i,:)./(r_3_mag(i,:)).^3); % in m/s^2
end
end
