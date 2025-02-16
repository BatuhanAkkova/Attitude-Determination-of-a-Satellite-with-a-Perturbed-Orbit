function s = Sun_Sensor_Modelling(jd,nu,omega,RAAN,inc,Period)
% Inputs:
% jd: Julian Date on Orbit
% nu: True Anomaly in rad
% omega: Argument of Perigee in deg
% RAAN: Right Asc. of Asc. Node in deg
% inc: Inclination in deg
% Period: Orbit's Period

% Output: Sun's Attitude from Reference Frame

s = zeros(Period,3); % Preallocation

for i = 1:Period

T_UTI = (jd(i)-2451545)/36525; % Number of Julian Centuries from epoch

% Mean Longitude
Longitude_Mean = 280.460 + 36000.771*T_UTI;
Longitude_Mean = deg2rad(Longitude_Mean);
Longitude_Mean = mod(Longitude_Mean,2*pi);

% Mean Anomaly
Mean_Anomaly = 357.5291092 + 35999.05034*T_UTI;
Mean_Anomaly = deg2rad(Mean_Anomaly);
Mean_Anomaly = mod(Mean_Anomaly,2*pi);

% Ecliptic Latitude
lambda_ecliptic = Longitude_Mean + 1.914666471*sin(Mean_Anomaly) + 0.019994643*sin(2*Mean_Anomaly);
lambda_ecliptic = deg2rad(lambda_ecliptic);
lambda_ecliptic = mod(mod(lambda_ecliptic,2*pi),pi/2);

% Obliquity of the ecliptic
E_oblq = 23.439291 - 0.0130042*T_UTI;
E_oblq = deg2rad(E_oblq);

% Attitude vector in geocentric equatorial frame (MOD)
s(i,:) = [cos(lambda_ecliptic) cos(E_oblq).*sin(lambda_ecliptic) sin(E_oblq).*sin(lambda_ecliptic)];

% DCM to convert s to spacecraft frame

if size(omega,1) > 1 % If variations are input:
theta = nu(i)-omega(i);
DCM = [cos(RAAN(i))*cos(theta)-sin(RAAN(i))*cos(inc(i))*sin(theta) -cos(RAAN(i))*sin(theta)-sin(RAAN(i))*cos(inc(i))*cos(theta) sin(RAAN(i))*sin(inc(i));...
    sin(RAAN(i))*cos(theta)+cos(RAAN(i))*cos(inc(i))*sin(theta) -sin(RAAN(i))*sin(theta)+cos(RAAN(i))*cos(inc(i))*cos(theta) -cos(RAAN(i))*sin(inc(i));...
    sin(inc(i))*sin(theta) sin(inc(i))*cos(theta) cos(inc(i))];
else
theta = nu(i)-omega;
DCM = [cos(RAAN)*cos(theta)-sin(RAAN)*cos(inc)*sin(theta) -cos(RAAN)*sin(theta)-sin(RAAN)*cos(inc)*cos(theta) sin(RAAN)*sin(inc);...
    sin(RAAN)*cos(theta)+cos(RAAN)*cos(inc)*sin(theta) -sin(RAAN)*sin(theta)+cos(RAAN)*cos(inc)*cos(theta) -cos(RAAN)*sin(inc);...
    sin(inc)*sin(theta) sin(inc)*cos(theta) cos(inc)];
end
s(i,:) = s(i,:)*DCM;
end