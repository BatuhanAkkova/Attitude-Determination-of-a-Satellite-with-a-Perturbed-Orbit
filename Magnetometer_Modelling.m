% Magnetometer Modelling
function B = Magnetometer_Modelling(LLA,jd,dt,rot_E,Period)

lat = LLA(:,1); % Latitude
lon = LLA(:,2); % Longitude
height = abs(LLA(:,3)); % Altitude

% Converting julian date to decimal year for igrfmagm function
Dec_Year = decyear(datetime(jd,'ConvertFrom','juliandate'));

% The Earth's Rotational Angle
%https://en.wikipedia.org/wiki/Sidereal_time#Earth_rotation_angle
ERA = 2*pi*(0.7790572732640 + 1.00273781191135448*(jd-2451545.0));
ERA = mod(ERA,2*pi);

gamma = rot_E.*dt + ERA; % angle to transform eci and ecef
gamma = mod(gamma,2*pi);

% Magnetic Field in NED Frame with IGRF-13
B_NED = igrfmagm(height,lat,lon,Dec_Year); % Magnetic field calculation with igrf-13

% Ned to Ecef transformation
[B_ECEF(:,1),B_ECEF(:,2),B_ECEF(:,3)] = ned2ecef(B_NED(:,1),B_NED(:,2),B_NED(:,3),lat,lon,height,wgs84Ellipsoid);

% Rotational Matrix From ECI to ECEF
% https://space.stackexchange.com/questions/38807/transform-eci-to-ecef
B = zeros(Period,3); % Preallocation for speed
for i=1:Period
    R_ECI2ECEF = [cos(gamma(i)) -sin(gamma(i)) 0; sin(gamma(i)) cos(gamma(i)) 0; 0 0 1];
    B(i,:) = B_ECEF(i,:)*R_ECI2ECEF'; % Magnetic Field in ECI
end
end