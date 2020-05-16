clear all;
close all;
clc;

%% Data
% Physical constants
Physics.R = 8.31446261815324;

% Vehicle data
Vehicle.mass = 1500;
Vehicle.J = 9.28951e+02;
Vehicle.S = 8.498096887876168;
Vehicle.d = 0.05;
Vehicle.CD = 1.5;
Vehicle.CL = 0.8;
Vehicle.CMalpha = -0.07;
Vehicle.CMq = -0.05;
Vehicle.CMdelta = 0.10;

Vehicle.Parachute.CD = 1.17;
Vehicle.Parachute.D = 40;
Vehicle.Parachute.S =  4 * pi * Vehicle.Parachute.D^2 / 8;
Vehicle.Parachute.mass = 0;

% Mars data
Mars.radius = 3397e3;
Mars.mu = 42830e9;

% Atmosphere data
Atm.rho0 = 0.02;
Atm.hs = 11100;
Atm.gamma = 1.29;
Atm.mean_molar_mass = 43.34e-3;
Atm.R = Physics.R / Atm.mean_molar_mass;

% Simulation parameters
params.parachute.at = 4000;

params.retro_rocket.at = 100; % meters
params.retro_rocket.force = 100; % Newton

params.r_des = params.parachute.at + Mars.radius;
params.v_des = 200;
params.rho_des = Atm.rho0*exp(-params.parachute.at/Atm.hs);
params.tau = 0.2;
params.zeta = 0.7;
params.wn = 20;
params.controlled = true;

%% Initial conditions
v_ini = 3.555201017579192e+03;
gamma_ini = -0.799040845502683;
h_ini = 120000;
phi_ini = deg2rad(0.0);
theta_ini = deg2rad(-80);
q_ini = deg2rad(0.0);

xi = [v_ini, gamma_ini, h_ini, phi_ini, theta_ini, q_ini];

%% Simulation
output = sim("model.slx", "RelTol", "1e-6", "StartTime", "0", "StopTime", "1000");

%% Assignation
t = output.v.time;
v = output.v.data;
gamma = output.gamma.data;
h = output.h.data;
phi = output.phi.data;
theta = output.theta.data;
q = output.q.data;

%% Intermediate calculations
r = h + Mars.radius;
g = Mars.mu / r.^2; % m/s^2 - Gravitationnal acceleration at r
rho =  Atm.rho0 * exp(-h/Atm.hs);
Pdyn = (1/2) .* rho .* v.^2; % Dynamic p ressure
alpha = theta - gamma;
M = MachNumber(v, h, Atm);

gamma_ref = zeros(size(t));

if (params.controlled)
    for i = 1:size(gamma_ref)
        if (h(i) > params.parachute_at)
            B = (Vehicle.S .* Vehicle.CD) ./ (Vehicle.mass);
            dv_aero = params.v_des - sqrt(v(i).^2+((2.*Mars.mu).*((1./params.r_des)-(1./r(i)))));
            gamma_ref(i) = asin((1/2).*B.*Atm.hs.*((params.rho_des - rho(i)) ./ log(1 + dv_aero./v(i))));
        else
            gamma_ref(i) = deg2rad(-90);
        end
    end
end

Daero = Pdyn .* Vehicle.S .* Vehicle.CD;
Laero = Pdyn .* Vehicle.S .* Vehicle.CL .* alpha;

%% Plotting
figure;
grid on;
plot(t, h);
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude function of time");

figure;
grid on;
plot(t, v);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
title("Velocity function of time");

figure;
grid on;
plot(h, v);
xlabel("Altitude (m)");
ylabel("Velocity (m/s)");
title("Velocity function of altitude");
set(gca,'Xdir','reverse');

figure;
grid on;
hold on;
plot(t, Laero);
plot(t, Daero);
legend("Laero", "Daero");
xlabel("Time (s)");
ylabel("Aerodynamic forces");
title("L_{aero} and D_{aero} function of time");

figure;
grid on;
hold on;
plot(t, rad2deg(gamma));
plot(t, rad2deg(gamma_ref));
legend("\gamma", "\gamma_{ref}");
title("Flight path angle \gamma function of time");
xlabel("Time (s)");
ylabel("Flight path angle (deg)");

figure;
grid on;
plot(h, rad2deg(gamma));
title("Flight path angle \gamma function of altitude");
xlabel("Altitude (m)");
ylabel("Flight path angle (deg)");
set(gca,'Xdir','reverse');

figure;
grid on;
hold on;
plot(t, theta);
plot(t, alpha);
legend("\theta", "\alpha");
title("\alpha and \theta function of time");
xlabel("Time (s)");
ylabel("Angles (deg)");

figure;
plot(t, M);
title("Mach number function of time");
xlabel("Time (s)");
ylabel("Mach number");
grid on;

figure;
plot(h, M);
title("Mach number function of altitude");
xlabel("Altitude");
ylabel("Mach number (m)");
set(gca,'Xdir','reverse');
grid on;