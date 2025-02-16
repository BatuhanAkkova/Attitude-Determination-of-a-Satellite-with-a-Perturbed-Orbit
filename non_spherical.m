% Non-Spherical Earth Perturbations
function agrav = non_spherical(Re,state,mu,Period)
% Inputs:
% Re: radius of Earth
% state: r vector
% mu: Gravitational parameter of Earth

% Output: acceleration due to non-spherical Earth in m/s^2

% Non-Spherical Earth Coefficients REF: EGM 2008
J2 = 1082.6267/10^6;
J3 = -2.5327/10^6;
J4 = -1.6196/10^6;
J5 = -0.2273/10^6;
J6 = 0.6083/10^8;

agrav = zeros(Period,3);

for i = 1:Period

r = sqrt(state(i,1)^2 + state(i,2)^2 + state(i,3)^2); % Magnitude of r

% Acceleration due to J2
a_J2 = [-3*J2*mu*Re^2*state(i,1)/(2*r^5)*(1 - 5*state(i,3)^2/r^2);...
    -3*J2*mu*Re^2*state(i,2)/(2*r^5)*(1 - 5*state(i,3)^2/r^2);...
    -3*J2*mu*Re^2*state(i,3)/(2*r^5)*(3 - 5*state(i,3)^2/r^2)];

% Acceleration due to J3
a_J3 = [-5*J3*mu*Re^3*state(i,1)/(2*r^7)*(3*state(i,3) - 7*state(i,3)^3/r^2);...
    -5*J3*mu*Re^3*state(i,2)/(2*r^7)*(3*state(i,3) - 7*state(i,3)^3/r^2);...
    -5*J3*mu*Re^3/(2*r^7)*(6*state(i,3)^2 - 7*state(i,3)^4/r^2 - 3*r^2/5)];

% Acceleration due to J4
a_J4 = [15*J4*mu*Re^4*state(i,1)/(8*r^7)*(1 - 14*(state(i,3))^2/r^2 + 21*(state(i,3))^4/r^4);...
    15*J4*mu*Re^4*state(i,2)/(8*r^7)*(1 - 14*(state(i,3))^2/r^2 + 21*(state(i,3))^4/r^4);...
    15*J4*mu*Re^4*state(i,3)/(8*r^7)*(5 - 70*(state(i,3))^2/(3*r^2) + 21*(state(i,3))^4/r^4)];

% Acceleration due to J5
a_J5 = [3*J5*mu*Re^5*state(i,1)*state(i,3)/(8*r^9)*(35 - 210*state(i,3)^2/r^2 + 231*state(i,3)^4/r^4);...
    3*J5*mu*Re^5*state(i,2)*state(i,3)/(8*r^9)*(35 - 210*state(i,3)^2/r^2 + 231*state(i,3)^4/r^4);...
    3*J5*mu*Re^5*state(i,3)*state(i,3)/(8*r^9)*(105 - 315*state(i,3)^2/r^2 + 231*state(i,3)^4/r^4)+15*J5*mu*Re^5/(8*r^7)];

% Acceleration due to J6
a_J6 = [-J6*mu*Re^6*state(i,1)/(16*r^9)*(35 - 945*state(i,3)^2/r^2 + 3465*state(i,3)^4/r^4 - 3003*state(i,3)^6/r^6);...
    -J6*mu*Re^6*state(i,2)/(16*r^9)*(35 - 945*state(i,3)^2/r^2 + 3465*state(i,3)^4/r^4 - 3003*state(i,3)^6/r^6);...
    -J6*mu*Re^6/(16*r^9)*(245 - 2205*state(i,3)^2/r^2 + 4851*state(i,3)^4/r^4 - 3003*state(i,3)^6/r^6)];

% Non-Spherical Acceleration
agrav(i,:) = a_J2 + a_J3 + a_J4 + a_J5 + a_J6;
end
end