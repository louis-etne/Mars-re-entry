function [Orbit_ECI] = orbites(a_m, e, M_deg, i_deg)

mu = 3.98618e14;
ID = 25544;
%% *Orbital Plane Coordinates*

n_rad_per_s = sqrt(mu/a_m^3);  % [rad/s] mean motion
n_deg_per_s = rad2deg(n_rad_per_s); % [deg/s] mean motion
M_rad = deg2rad(M_deg);
E_rad = M_rad; 
dE = 99999;
eps = 1e-8; % [rad] control precision of Newton's method solution
while (abs(dE) > eps)
    dE = (E_rad - e * sin(E_rad) - M_rad)/(1 - e * cos(E_rad));
    E_rad = E_rad -  dE;
end

E_deg_epoch = rad2deg(E_rad); 
Omega_deg = 0;
%% *Rotate To ECI*
Rz_Omega = [ ...
    [cosd(Omega_deg) sind(Omega_deg) 0]; ...
    [-sind(Omega_deg) cosd(Omega_deg) 0]; ...
    [0 0 1]];
Rx_i = [ ...
    [1 0 0]; ...
    [0 cosd(i_deg) sind(i_deg)]; ...
    [0 -sind(i_deg) cosd(i_deg)]];
Rz_omega = [ ...
    [cosd(Omega_deg) sind(Omega_deg) 0]; ...
    [-sind(Omega_deg) cosd(Omega_deg) 0]; ...
    [0 0 1]];


Evals = 0:1:360.0; % [deg] values of the eccentric anomaly around orbit 
Orbit_p = a_m*(cosd(Evals)-e); % [m] orbit positions
Orbit_q = a_m*sqrt(1 - e^2)*sind(Evals); % [m] orbit positions
deltaT_s = ((Evals-E_deg_epoch) - e*sind(Evals-E_deg_epoch))/n_deg_per_s; % [s] time since epoch along orbit
Orbit_ECI = zeros(numel(deltaT_s),3);
for ipt = 1:size(Orbit_ECI,1)
    r_pq = [Orbit_p(ipt) Orbit_q(ipt) 0]';
    Orbit_ECI(ipt,:) = [inv(Rz_Omega)*inv(Rx_i)*inv(Rz_omega)*r_pq]'; %[Rz_Omega*Rx_i*Rz_omega*r_pq]';
end
end

