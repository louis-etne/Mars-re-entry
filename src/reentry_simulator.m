clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat');

coefs.CA = load('fit_areo_axial_coef.mat');
coefs.CN = load('fit_areo_norm_coef.mat');

coeff_orbit_entry = load('orbit.mat');


%% Initial conditions
v_ini = coeff_orbit_entry.V_entry; %  m/s - Initial velocity at h_ini
gamma_ini = -coeff_orbit_entry.flight_path_angle; % rad - Angle between horizontal plane and the velocity vector - Initial flight path angle
h_ini = 60000; % m - Initial altitude
phi_ini = deg2rad(0.0); % rad - Downrange
theta_ini = deg2rad(-80); % rad - Angle between horizontal and the capsule axis
q_ini = deg2rad(0.0); % rad/s - Angular speed between horizontal and the capsule axis 

xi = [v_ini, gamma_ini, h_ini, phi_ini, theta_ini, q_ini];

%% Simulation parameters
sim_params.parachute_at = 8000; % m - Open parachute at 8km
sim_params.parachute_mach = 1.8; % Max mach number for parachute
sim_params.r_des = sim_params.parachute_at + Mars.radius * 1000;
sim_params.v_des = sim_params.parachute_mach * SpeedOfSound(sim_params.parachute_at, Atm);
sim_params.rho_des = Density(sim_params.parachute_at, Atm);

% Controller parameters
sim_params.tau = 0.2;
sim_params.zeta = 0.7;
sim_params.wn = 20;

sim_params.controlled = false;
ts = [0 300]; % s - Time span, [t0 tf]

%% Simulation
options = odeset('reltol', 1e-8);
tic
[t, x] = ode45(@(t,y)Dynamics(t, y, Mars, Atm, Vehicle, coefs, sim_params), ts, xi, options);
toc
x = real(x);

%% Assignations
v = x(:, 1);
gamma = x(:, 2);
h = x(:, 3);
phi = x(:, 4);
theta = x(:, 5);
q = x(:, 6);

save('Data/flight.mat', 't', 'x', 'v', 'gamma', 'h', 'phi', 'theta', 'q');

%% Simulation Computations
alpha = theta - gamma;
M = MachNumber(v, h, Atm);
rho = Density(h, Atm); % kg/m^3 - Air density at h
Pdyn = (1/2) .* rho .* v.^2; % Dynamic pressure

CA = [];
CN = [];

for i = 1:size(M)
    CA = [CA; AxialForceCoef(M(i), alpha(i), coefs.CA)];
    CN = [CN; NormalForceCoef(M(i), alpha(i), coefs.CN)];
end

CL = CN;
CD = CA;

Daero = Pdyn .* Vehicle.S .* CD;
Laero = Pdyn .* Vehicle.S .* CL .* alpha;

%% Plots
figure;
hold on;
plot(t, h);
grid on;
title('Altitude vs time');
xlabel('Time (s)');
ylabel('Altitude (m)');

figure;
hold on;
plot(t, v);
grid on;
title('Velocity vs time');
xlabel('Time (s)');
ylabel('Velocity (m.s^{-1})');

figure;
ax = axes;
hold on;
plot(h, v);
grid on;
title('Velocity versus altitude');
xlabel('Altitude (m)');
ylabel('Velocity (m.s^{-1})');
set(ax, 'XDir', 'reverse');

figure;
hold on;
grid on;
plot(t, CL);
plot(t, CD);
legend('Lift coefficient', 'Drag coefficient');
title('Aerodynamic coefficients versus time');
xlabel('Time (s)');
ylabel('Aerodynamic coefficients');

figure;
hold on;
grid on;
plot(t, Laero);
plot(t, Daero);
legend('Lift', 'Drag');
title('D_{aero} and L_{aero} versus time');
xlabel('Time (s)');
ylabel('Aerodynamic forces');

figure;
hold on;
grid on;
plot(t, rad2deg(alpha));
plot(t, rad2deg(theta));
legend('\alpha', '\theta');
title('\alpha and \theta versus time');
xlabel('Time (s)');
ylabel('Angles (deg)');
