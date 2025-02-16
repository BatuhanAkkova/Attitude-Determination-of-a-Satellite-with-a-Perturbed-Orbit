initial_values;

        % TIME SCRIPT

[utc,jd] = Time(Period);

        % KEPLER'S EQUATION SCRIPT

E = kepler_eq(e,Me,Period); % Eccentric Anomaly Calc.

        % STATE VECTOR SCRIPT

[r.kepler,v.kepler,nu] = Kepler_Elements_to_State(e,inc,RAAN,omega,E,Period,mu,h);
LLA.kepler = eci2lla(r.kepler,utc); % Latitude, Longitude, Altitude in geodetic

        % CALCULATING KINEMATIC AND DYNAMIC EQUATIONS

q_true = EOM(I,L,w0,q0,Period);

% Corresponding direction cosine matrix:
A_true = zeros(Period*3,3); % Preallocation of True Rotational Matrix
ii = 1;
for i=1:3:Period*3
A_true(i:i+2,:) = quat2dcm(q_true(ii,:));
ii = ii+1;
end
clear i ii

        % ATMOSPHERIC DENSITY (via Exponential Method)
% rho0: 3.019e-15 km/m3 Nominal Density
% h0: 1000e3 m Base Altitude
% H: 268e3 m Scale Height

rho = 3.019E-15*exp(-(LLA.kepler(:,3)-1E6)/268E3); 

        % PERTURBATIONS

a.grav = non_spherical(Re,r.kepler,mu,Period);
a.drag = atm_drag(rho,v.kepler,m,A);
a.third = Third_Body(r.kepler,jd,Period);
a.solar = solar_p(m,A,r.kepler,jd);

a.total = a.grav + a.drag + a.third + a.solar;

% Computing State Vector with Encke using various accelerations
% Total Pert.
SOLN.total = Encke(r.kepler,v.kepler,mu,h,Period,a.total);
r.total = SOLN.total(2:end,1:3);
v.total = SOLN.total(2:end,4:6);
LLA.total = eci2lla(r.total,utc);

% Only Grav.
SOLN.grav = Encke(r.kepler,v.kepler,mu,h,Period,a.grav);
r.grav = SOLN.grav(2:end,1:3);
v.grav = SOLN.grav(2:end,4:6);
LLA.grav = eci2lla(r.grav,utc);

% Only Drag
SOLN.drag = Encke(r.kepler,v.kepler,mu,h,Period,a.drag);
r.drag = SOLN.drag(2:end,1:3);
v.drag = SOLN.drag(2:end,4:6);
LLA.drag = eci2lla(r.drag,utc);

% Only Third Body
SOLN.third = Encke(r.kepler,v.kepler,mu,h,Period,a.third);
r.third = SOLN.third(2:end,1:3);
v.third = SOLN.third(2:end,4:6);
LLA.third = eci2lla(r.third,utc);

% Only Solar
SOLN.solar = Encke(r.kepler,v.kepler,mu,h,Period,a.solar);
r.solar = SOLN.solar(2:end,1:3);
v.solar = SOLN.solar(2:end,4:6);
LLA.solar = eci2lla(r.solar,utc);

% Only the Max.
SOLN.max = Encke(r.kepler,v.kepler,mu,h,Period,a.grav);
r.max = SOLN.max(2:end,1:3);
v.max = SOLN.max(2:end,4:6);
LLA.max = eci2lla(r.max,utc);

        % Calculating Variation of Parameters

% Total Pert.
[e_pert.total,inc_pert.total,RAAN_pert.total,omega_pert.total,h_pert.total,nu_pert.total] = State_to_Kepler_Elements(r.total,v.total,mu,Period);

% Only Grav.
[e_pert.grav,inc_pert.grav,RAAN_pert.grav,omega_pert.grav,h_pert.grav,nu_pert.grav] = State_to_Kepler_Elements(r.grav,v.grav,mu,Period);

% Only Drag
[e_pert.drag,inc_pert.drag,RAAN_pert.drag,omega_pert.drag,h_pert.drag,nu_pert.drag] = State_to_Kepler_Elements(r.drag,v.drag,mu,Period);

% Only Third Body
[e_pert.third,inc_pert.third,RAAN_pert.third,omega_pert.third,h_pert.third,nu_pert.third] = State_to_Kepler_Elements(r.third,v.third,mu,Period);

% Only Solar
[e_pert.solar,inc_pert.solar,RAAN_pert.solar,omega_pert.solar,h_pert.solar,nu_pert.solar] = State_to_Kepler_Elements(r.solar,v.solar,mu,Period);

% Only the Max.
[e_pert.max,inc_pert.max,RAAN_pert.max,omega_pert.max,h_pert.max,nu_pert.max] = State_to_Kepler_Elements(r.max,v.max,mu,Period);

        % CALCULATING SUN SENSOR MODEL

s.total = Sun_Sensor_Modelling(jd,nu_pert.total,omega_pert.total,RAAN_pert.total,inc_pert.total,Period); % For Sensor Reference

% For QUEST Input, change between various perturbations:
s.grav = Sun_Sensor_Modelling(jd,nu_pert.grav,omega_pert.grav,RAAN_pert.grav,inc_pert.grav,Period); % Only gravity
s.drag = Sun_Sensor_Modelling(jd,nu_pert.drag,omega_pert.drag,RAAN_pert.drag,inc_pert.drag,Period); % Only drag
s.solar = Sun_Sensor_Modelling(jd,nu_pert.solar,omega_pert.solar,RAAN_pert.solar,inc_pert.solar,Period); % Only solar
s.third = Sun_Sensor_Modelling(jd,nu_pert.third,omega_pert.third,RAAN_pert.third,inc_pert.third,Period); % Only third
s.max = Sun_Sensor_Modelling(jd,nu_pert.max,omega_pert.max,RAAN_pert.max,inc_pert.max,Period); % Only the max.
s.kepler = Sun_Sensor_Modelling(jd,nu,deg2rad(omega),deg2rad(RAAN),deg2rad(inc),Period); % No pert.

% Reference*A_true = Body Frame
s_meas = zeros(Period,3); % Preallocation
ii = 1;
for i=1:Period
s_meas(i,:) = s.total(i,:)*A_true(ii:ii+2,:) + randn*deviation_sun; % Measured Sun Sensor Attitude in Body Frame
ii = ii+3;
end
clear i ii

        % CALCULATING MAGNETOMETER MODEL        
B.total = Magnetometer_Modelling(LLA.total,jd,dt,rot_E,Period); % For Sensor Reference

% For QUEST Input, change between various perturbations:
B.grav = Magnetometer_Modelling(LLA.grav,jd,dt,rot_E,Period); % Only gravity
B.drag = Magnetometer_Modelling(LLA.drag,jd,dt,rot_E,Period); % Only drag
B.solar = Magnetometer_Modelling(LLA.solar,jd,dt,rot_E,Period); % Only solar
B.third = Magnetometer_Modelling(LLA.third,jd,dt,rot_E,Period); % Only third body
B.max = Magnetometer_Modelling(LLA.max,jd,dt,rot_E,Period); % Only max
B.kepler = Magnetometer_Modelling(LLA.kepler,jd,dt,rot_E,Period); % No pert.

% Reference*A_true = Body Frame
B_meas = zeros(Period,3); % Preallocation
ii = 1;
for i=1:Period
B_meas(i,:) = B.total(i,:)*A_true(ii:ii+2,:) + randn*deviation_mag; % Measured Magnetic Field in body frame
ii = ii+3;

% Unit Vector:

B.total(i,:) = B.total(i,:) / norm(B.total(i,:));
B.grav(i,:) = B.grav(i,:) / norm(B.grav(i,:));
B.drag(i,:) = B.drag(i,:) / norm(B.drag(i,:));
B.solar(i,:) = B.solar(i,:) / norm(B.solar(i,:));
B.third(i,:) = B.third(i,:) / norm(B.third(i,:));
B.max(i,:) = B.max(i,:) / norm(B.max(i,:));
B.kepler(i,:) = B.kepler(i,:) / norm(B.kepler(i,:));

B_meas(i,:) = B_meas(i,:) / norm(B_meas(i,:));
end
clear i ii

        % CALCULATING ATTITUDE MATRIX VIA QUEST METHOD

% TOTAL - All the Perturbations Included
[A_opt.total,q_opt.total] = QUEST(B.total,B_meas,s.total,s_meas,Period);

% Multiplying (-) if needed (for quaternion errors)
for i=1:Period
for ii=1:4
if q_true(i,ii) < 0
    q_opt.total(i,ii) = -abs(q_opt.total(i,ii));
elseif q_true(i,ii) > 0
    q_opt.total(i,ii) = abs(q_opt.total(i,ii));
end
end
end
clear i ii

% GRAV - Only the Gravitational Perturbations (e.g. J2) Included)
[A_opt.grav,q_opt.grav] = QUEST(B.grav,B_meas,s.grav,s_meas,Period);

% Multiplying (-) if needed (for quaternion errors)
for i=1:Period
for ii=1:4
if q_true(i,ii) < 0
    q_opt.grav(i,ii) = -abs(q_opt.grav(i,ii));
elseif q_true(i,ii) > 0
    q_opt.grav(i,ii) = abs(q_opt.grav(i,ii));
end
end
end
clear i ii

% DRAG - Only the Atmospheric Drag Perturbation Included
[A_opt.drag,q_opt.drag] = QUEST(B.drag,B_meas,s.drag,s_meas,Period);

% Multiplying (-) if needed (for quaternion errors)
for i=1:Period
for ii=1:4
if q_true(i,ii) < 0
    q_opt.drag(i,ii) = -abs(q_opt.drag(i,ii));
elseif q_true(i,ii) > 0
    q_opt.drag(i,ii) = abs(q_opt.drag(i,ii));
end
end
end
clear i ii

% SOLAR - Only the Solar Radiation Pressure Included
[A_opt.solar,q_opt.solar] = QUEST(B.solar,B_meas,s.solar,s_meas,Period);

% Multiplying (-) if needed (for quaternion errors)
for i=1:Period
for ii=1:4
if q_true(i,ii) < 0
    q_opt.solar(i,ii) = -abs(q_opt.solar(i,ii));
elseif q_true(i,ii) > 0
    q_opt.solar(i,ii) = abs(q_opt.solar(i,ii));
end
end
end
clear i ii

% THIRD - Only the Third Body (the Moon) Perturbation Included
[A_opt.third,q_opt.third] = QUEST(B.third,B_meas,s.third,s_meas,Period);

% Multiplying (-) if needed (for quaternion errors)
for i=1:Period
for ii=1:4
if q_true(i,ii) < 0
    q_opt.third(i,ii) = -abs(q_opt.third(i,ii));
elseif q_true(i,ii) > 0
    q_opt.third(i,ii) = abs(q_opt.third(i,ii));
end
end
end
clear i ii

%  MAX - Only the Maximum Effective (probably Grav.) Perturbation Included
[A_opt.max,q_opt.max] = QUEST(B.max,B_meas,s.max,s_meas,Period);

% Multiplying (-) if needed (for quaternion errors)
for i=1:Period
for ii=1:4
if q_true(i,ii) < 0
    q_opt.max(i,ii) = -abs(q_opt.max(i,ii));
elseif q_true(i,ii) > 0
    q_opt.max(i,ii) = abs(q_opt.max(i,ii));
end
end
end
clear i ii

% NO PERT. - No Perturbations are Included
[A_opt.kepler,q_opt.kepler] = QUEST(B.kepler,B_meas,s.kepler,s_meas,Period);

% Multiplying (-) if needed (for quaternion errors)
for i=1:Period
for ii=1:4
if q_true(i,ii) < 0
    q_opt.kepler(i,ii) = -abs(q_opt.kepler(i,ii));
elseif q_true(i,ii) > 0
    q_opt.kepler(i,ii) = abs(q_opt.kepler(i,ii));
end
end
end
clear i ii

angle_quest = zeros(1,Period); % Preallocation

% The angle between B and s for the possible errors 
for i = 1:size(B.total,1)
    angle_quest(i) = acosd(dot(s.total(i,:), B.total(i,:)));
end

% Change in quaternions
delta_q = q_opt.total - q_true;

% PLOTTING THE OUTPUTS
N=Period;

% PLOTTING THE KEPLER ORBIT
% Earth as a sphere
[X,Y,Z] = sphere; 
X = X*Re;
Y = Y*Re;
Z = Z*Re;

figure(1);
plot3(r.kepler(1:Period-1,1),r.kepler(1:Period-1,2),r.kepler(1:Period-1,3),'black','LineWidth',2) % Orbit
title('Kepler Orbit')
xlabel('X-axis (m)')
ylabel('Y-axis (m)')
zlabel('Z-axis (m)')
grid on
hold on
surf(X,Y,Z) % Earth
colormap([0 1 1])
hold off

% PLOTTING THE PERTURBED ORBIT
% Earth as a sphere
[X,Y,Z] = sphere; 
X = X*Re;
Y = Y*Re;
Z = Z*Re;

figure(2);
plot3(r.total(:,1),r.total(:,2),r.total(:,3),'black','LineWidth',2) % Orbit
title('Perturbed Orbit')
xlabel('X-axis (m)')
ylabel('Y-axis (m)')
zlabel('Z-axis (m)')
grid on
hold on
surf(X,Y,Z) % Earth
colormap([0 1 1])
hold off

% PLOTTING TRUE QUATERNION
figure(3);
t5 = tiledlayout(2,2,"TileSpacing","compact");
title(t5,'True Quaternion')
xlabel(t5,'Time (s)')
ylabel(t5,'Attitude')

nexttile
plot(dt,q_true(:,1))
title('q0')
grid on

nexttile
plot(dt,q_true(:,2))
title('q1')
grid on

nexttile
plot(dt,q_true(:,3))
title('q2')
grid on

nexttile
plot(dt,q_true(:,4))
title('q3')
grid on

% PLOTTING TRUE QUATERNION AND ESTIMATED QUATERNION
figure(4);
t5 = tiledlayout(2,2,"TileSpacing","compact"); 
title(t5,'q-true vs q-opt')
xlabel(t5,'Time (s)')
ylabel(t5,'Attitude Value')

nexttile
plot(dt,q_true(1:N,1))
hold on
plot(dt,q_opt.total(:,1))
hold off
legend('q-true','q-opt')
title('q0')
grid on

nexttile
plot(dt,q_true(1:N,2))
hold on
plot(dt,q_opt.total(1:N,2))
hold off
legend('q-true','q-opt')
title('q1')
grid on

nexttile
plot(dt,q_true(1:N,3))
hold on
plot(dt,q_opt.total(1:N,3))
hold off
legend('q-true','q-opt')
title('q2')
grid on

nexttile
plot(dt,q_true(1:N,4))
hold on
plot(dt,q_opt.total(1:N,4))
hold off
legend('q-true','q-opt')
title('q3')
grid on

% PLOTTING PERTURBATIONS ONE-BY-ONE

figure(5);
plot(dt,a.grav)
title('Non-Spherical Earth')
xlabel('Time (s)')
ylabel('Acceleration (m/s2)')
grid on

figure(6);
plot(dt,a.drag)
title('Atmospheric Drag')
xlabel('Time (s)')
ylabel('Acceleration (m/s2)')
grid on

figure(7);
plot(dt,a.third)
title('Third-Body G. Effect')
xlabel('Time (s)')
ylabel('Acceleration (m/s2)')
grid on

figure(8);
plot(dt,a.solar)
title('Solar Radiation Pressure')
xlabel('Time (s)')
ylabel('Acceleration (m/s2)')
grid on

% PLOTTING VARIATION OF PARAMETERS

figure(9);
t7 = tiledlayout(3,2,"TileSpacing","compact");
title(t7,'Variation of Parameters')
xlabel(t7,'Time')

nexttile
plot(dt,vecnorm(r.kepler,2,2)-Re)
title('Altitude')
grid on

nexttile
plot(dt,e_pert.total)
title('Eccentricity')
grid on

nexttile
plot(dt,h_pert.total)
title('Angular Momentum')
grid on

nexttile
plot(dt,inc_pert.total)
title('Inclination')
grid on

nexttile
plot(dt,RAAN_pert.total)
title('Right Asc.')
grid on

nexttile
plot(dt,omega_pert.total)
title('Arg. of Perigee')
grid on

% PLOTTING TOTAL PERTURBATION ACCELERATION
figure(10);
plot(dt,a.total)
title('Total Perturbation Acceleration')
xlabel('Time (s)')
ylabel('Acceleration (m/s2)')
grid on

% PLOTTING DELTA_Q AND ANGLE BETWEEN B AND s
figure(11);
t11 = tiledlayout(3,2);
title(t11,'Delta-q and Angle Between Measurement Vectors')
xlabel(t11,'Time (s)')

nexttile
plot(dt,angle_quest)
grid on

nexttile
plot(dt,delta_q(:,1))
grid on

nexttile
plot(dt,delta_q(:,2))
grid on

nexttile
plot(dt,delta_q(:,3))
grid on

nexttile
plot(dt,delta_q(:,4))
grid on

% RMSE TABLE

err.total = q_opt.total - q_true;
err.grav = q_opt.grav - q_true;
err.drag = q_opt.drag - q_true;
err.third = q_opt.third - q_true;
err.solar = q_opt.solar - q_true;
err.max = q_opt.max - q_true;
err.kepler = q_opt.kepler - q_true;

Table = [rms(err.total); rms(err.grav); rms(err.drag); rms(err.third); rms(err.solar); rms(err.max); rms(err.kepler)];

