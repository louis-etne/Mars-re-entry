clear all;
close all;
clc;

%% Data
% Physical constants
Physics.R = 8.31446261815324;

% Vehicle data
Vehicle.mass = 1443;
Vehicle.J = 9.28951e+02;
Vehicle.S = 8.498096887876168;
Vehicle.d = 0.05;
Vehicle.CD = 1.5;
Vehicle.CL = 0.8;
Vehicle.CMalpha = -0.07;
Vehicle.CMq = -0.05;
Vehicle.CMdelta = 0.10;
Vehicle.r_n = 1.59; % meters - nose radius

Vehicle.Parachute.CD = 1.75;
Vehicle.Parachute.D = 40;
Vehicle.Parachute.S =  4 * pi * Vehicle.Parachute.D^2 / 8;

% Additional masses
Vehicle.Parachute.mass = 62;
Vehicle.propellant = 80;
Vehicle.Heatshield.surface = 23.226; % m^2
% Vehicle.Heatshield.wetted_surface = 9.7; % m^2
% calculated after running a first simulation. 14.8201 is the specific
% mass of the heatshield
Vehicle.Heatshield.mass = 369;
                                                            
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

params.retro_rockets.at = 100; % m
params.retro_rockets.max_thrust = 1500*3; % N
params.retro_rockets.exhaust_velocity = 3000;   % m/s
params.retro_rockets.propellant = Vehicle.propellant;   % kg
params.retro_rockets.max_flow_rate = params.retro_rockets.max_thrust /...
                                params.retro_rockets.exhaust_velocity  % kg/s
params.retro_rockets.min_burn_duration = params.retro_rockets.propellant /...
                                    params.retro_rockets.max_flow_rate % s
                            
params.r_des = params.parachute.at + Mars.radius;
params.v_des = 200;
params.rho_des = Atm.rho0*exp(-params.parachute.at/Atm.hs);
params.tau = 0.2;
params.zeta = 0.7;
params.wn = 20;
params.controlled = true;

%% Initial conditions
v_ini = 3.555201017579192e+03;
% gamma_ini = -0.799040845502683; % wrong gamma 
gamma_ini = -0.058744571344473;
h_ini = 120000;
phi_ini = deg2rad(0.0);
theta_ini = deg2rad(-80);
q_ini = deg2rad(0.0);
m_ini = Vehicle.mass;
consumed_prop_ini = 0;
xi = [v_ini, gamma_ini, h_ini, phi_ini, theta_ini, q_ini, consumed_prop_ini, params.retro_rockets.max_thrust];

%% Simulation
output = sim("model.slx", "RelTol", "1e-6", "StartTime", "0");

%% Assignation
t = output.v.time;
v = output.v.data;
gamma = output.gamma.data;
h = output.h.data;
phi = output.phi.data;
theta = output.theta.data;
q = output.q.data;
mass = output.mass.data;
mass_capsule = output.mass_capsule.data;
heatflux = output.heatflux.data;
integrated_heatflux = output.integrated_heatflux.data;
acc = output.acc.data;
thrust = output.thrust.data;
% consumed_propellant = output.consumed_propellant.data;

%% Heatshield mass

peak_heatflux = max(integrated_heatflux);
heatshield_specific_mass = 0.24*peak_heatflux+0.29*sqrt(peak_heatflux)+11.3;
mass_heatshield = heatshield_specific_mass * Vehicle.Heatshield.surface;
ang_phi = 0:0.1:40;
q_c = peak_heatflux.*cosd(ang_phi).^1.2;
mass_h = 0;
for i=q_c
    mass_h = mass_h + 0.1*(0.24*i+0.29*sqrt(i)+11.3);
end
fun = @(x) (0.24 * (peak_heatflux.*cosd(x).^1.2) + 0.29 * sqrt(peak_heatflux.*cosd(x).^1.2)).*sin(x) + 11.3;

m_tot = Vehicle.r_n * 2 * pi * integral(fun, 0, 20);

%% Intermediate calculations

r = h + Mars.radius;
g = Mars.mu ./ r.^2; % m/s^2 - Gravitationnal acceleration at r
rho =  Atm.rho0 * exp(-h/Atm.hs);
Pdyn = (1/2) .* rho .* v.^2; % Dynamic p ressure
alpha = theta - gamma;
M = MachNumber(v, h, Atm);

gamma_ref = zeros(size(t));

if (params.controlled)
    for i = 1:size(gamma_ref)
        if (h(i) > params.parachute.at)
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
temp = (Vehicle.S .* (1000 * 3600 .* integrated_heatflux))./(mass_heatshield .* 750);

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
subplot(2, 3, 1)
plot(h, v);

xlabel("Altitude (m)");
ylabel("Velocity (m/s)");
title("Velocity function of altitude");
set(gca,'Xdir','reverse');
grid on;

figure;
grid on;
hold on;
plot(t(100: length(t)), Laero(100: length(t)));
plot(t(100: length(t)), Daero(100: length(t)));
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
plot(t, rad2deg(theta));
plot(t, rad2deg(alpha));
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

figure;
plot(mass_capsule, h());
title("Capsule mass function of altitude", 'Fontsize', 14);
xlabel("Mass (kg)", 'Fontsize', 14);
ylabel("Altitude (m)", 'Fontsize', 14);
set(gca,'Xdir','reverse');
grid on;

figure
plot(mass_capsule(1528:length(h)), h(1528:length(h)))
title("Mass variation during lander descent", 'Fontsize', 14);
xlabel("Mass (kg)", 'Fontsize', 14);
ylabel("Altitude (m)", 'Fontsize', 14);
set(gca,'Xdir','reverse');
grid on;

figure
plot(t(1531:length(h)), mass_capsule(1531:length(h)))
title("Mass variation during lander descent", 'Fontsize', 14);
ylabel("Mass (kg)", 'Fontsize', 14);
xlabel("Time(s)", 'Fontsize', 14);
grid on;


figure;
subplot(2, 1, 1);
plot(t, heatflux);
title("Heatflux function of time")
xlabel("Time (s)");
ylabel("Heatflux (kW/m^2)")
grid on;
subplot(2, 1, 2);
plot(t, integrated_heatflux, t, temp);
title("Integrated heatflux function of time")
xlabel("Time (s)");
ylabel("Heatflux (kWh/m^2)")
grid on;

figure;
plot(t, temp);
title("Temperature of heat shield function of time")
xlabel("Time (s)");
ylabel("Temperature (K)");
grid on;


figure;
% subplot(2, 1, 1)
plot(h(1531:length(h)), thrust(1531:length(h)), h(1531:length(h)), ...
    mass_capsule(1531:length(h)).*g(1531:length(h)));
title("Thrust function of altitude", 'Fontsize', 14)
xlabel("Altitude (m)", 'Fontsize', 14);
ylabel("Force (N)", 'Fontsize', 14);
legend("Thrust", "m\timesg");
set(gca,'Xdir','reverse');
grid on;
% subplot(2, 1, 2)
% plot(h(1528:length(h)), acc(1528:length(h)));
% title("Acceleration function of altitude", 'Fontsize', 14)
% xlabel("Altitude (m)", 'Fontsize', 14);
% ylabel("Thrust (N)", 'Fontsize', 14);
% set(gca,'Xdir','reverse');
% grid on;


% AnimatedLine
% lat = zeros(length(t));
% lon = linspace(1, 2, length(t));
% figure;
% grid on;
% ax = gca;
% % axis([-1, 1, 1, 2, 0, 120e3]);
% xlabel("Latitude (°)");
% ylabel("Longitude (°)");
% zlabel("Altitude (m)");
% ax.Clipping = 'off';
% title("Altitude function of time");
% curve = animatedline('Color', 'r');
% view(45, 45)
% for i = 1:length(t)
%    addpoints(curve, lat(i), lon(i), h(i));
%    drawnow;
% end

%%
alt = zeros(int16(0.1*length(h)), 1);
downrange = zeros(length(alt), 1);
time = zeros(length(alt), 1);
i = 1;
for j = 1:length(alt)
    alt(j) = h(i);
    downrange(j) = phi(i);
    time(j) = t(i);
    if i > 21139
        i = i + 1;
    else
        i = i + 10;
    end
end
% alt = h;
save('../src/Data/altitude.mat', 'alt', 'downrange', 'time');
