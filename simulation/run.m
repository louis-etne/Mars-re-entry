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
index_parachute = find(h <= params.parachute.at, 1);
t_parachute = floor(t(index_parachute));
vel_parachute = floor(v(index_parachute));

index_rockets = find(h <= params.retro_rockets.at, 1);
t_rockets = floor(t(index_rockets));

index_landing = length(h);
t_landing = floor(t(index_landing));

%% Acceleration function of time
figure
plot(t, abs(acc));
xlabel("Time (s)");
ylabel("Acceleration (m/s^2)")
title("Acceleration function of time");
grid on;
axis([0 t_parachute 0 100]);
%% Altitude function of time
figure;
subplot(2, 3, [1, 3]);
plot(t, h);
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude function of time");
grid on;
axis([0 t_landing 0 120000]);

% subplot(2, 3, 4);
%%
figure
plot(t, h);
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude function of time while gliding");
axis([0 t_parachute 4000 120000]);
grid on;
%%
subplot(2, 3, 5);
plot(t, h);
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude function of time with parachute");
axis([t_parachute t_rockets 150 4000]);
grid on;

subplot(2, 3, 6);
plot(t, h);
xlabel("Time (s)");
ylabel("Altitude (m)");
title("Altitude function of time with retrorockets");
axis([t_rockets t_landing 0 150]);
grid on;

%% Velocity function of time
figure;
subplot(2, 3, [1, 3]);
plot(t, v);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
title("Velocity function of time");
grid on;
axis([0 t_landing 0 3700]);

subplot(2, 3, 4);
plot(t, v);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
title("Velocity function of time while gliding");
axis([0 t_parachute vel_parachute 3700]);
grid on;

subplot(2, 3, 5);
plot(t, v);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
title("Velocity function of time with parachute");
axis([t_parachute t_rockets 10 226.2]);
grid on;

subplot(2, 3, 6);
plot(t, v);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
title("Velocity function of time with retrorockets");
axis([t_rockets t_landing 0 10]);
grid on;

%% Velocity function of altitude
figure;
subplot(2, 3, [1, 3]);
plot(h, v);
xlabel("Altitude (m)");
ylabel("Velocity (m/s)");
title("Velocity function of altitude");
axis([0 120000 0 3700]);
set(gca,'Xdir','reverse');
grid on;

subplot(2, 3, 4);
plot(h, v);
xlabel("Altitude (m)");
ylabel("Velocity (m/s)");
title("Velocity function of altitude while gliding");
axis([params.parachute.at 120000 vel_parachute 3700]);
set(gca,'Xdir','reverse');
grid on;

subplot(2, 3, 5);
plot(h, v);
xlabel("Altitude (m)");
ylabel("Velocity (m/s)");
title("Velocity function of altitude with parachute");
axis([params.retro_rockets.at 4000 0 vel_parachute]);
set(gca,'Xdir','reverse');
grid on;

subplot(2, 3, 6);
plot(h, v);
xlabel("Altitude (m)");
ylabel("Velocity (m/s)");
title("Velocity function of altitude with retrorockets");
axis([0 params.retro_rockets.at 0 10]);
set(gca,'Xdir','reverse');
grid on;

%% Laero Daero function of time
figure;
subplot(2, 2, [1, 2]);
hold on;
plot(t, Laero);
plot(t, Daero);
legend("L", "D");
xlabel("Time (s)");
ylabel("Aerodynamic forces");
title("L and D function of time");
grid on;
axis([0 t_landing -100 80000]);

subplot(2, 2, 3);
hold on;
plot(t, Laero);
plot(t, Daero);
legend("L", "D");
xlabel("Time (s)");
ylabel("Aerodynamic forces");
title("L and D function of time while gliding");
grid on;
axis([0 t_parachute -100 80000]);

subplot(2, 2, 4);
hold on;
plot(t, Laero);
plot(t, Daero);
legend("L", "D");
xlabel("Time (s)");
ylabel("Aerodynamic forces");
title("L and D function of time the rest of the descent");
grid on;
axis([t_parachute t_landing 0 20]);

%% L and D function of altitude
figure;
subplot(2, 2, [1, 2]);
hold on;
plot(h, Laero);
plot(h, Daero);
legend("L", "D", 'Location', 'northwest');
xlabel("Altitude (m)");
ylabel("Aerodynamic forces");
title("L and D function of altitude");
set(gca,'Xdir','reverse');
grid on;
axis([0 120000 0 190000]);

subplot(2, 2, 3);
hold on;
plot(h, Laero);
plot(h, Daero);
legend("L", "D", 'Location','northwest');
xlabel("Altitude (m)");
ylabel("Aerodynamic forces");
title("L and D function of altitude while gliding");
set(gca,'Xdir','reverse');
grid on;
axis([4000 120000 0 80000]);

subplot(2, 2, 4);
hold on;
plot(h, Laero);
plot(h, Daero);
legend("L", "D");
xlabel("Altitude (m)");
ylabel("Aerodynamic forces");
title("L and D function of altitude the rest of the descent");
set(gca,'Xdir','reverse');
grid on;
axis([0 4000 0 20]);

%% L and D function of velocity
figure;
subplot(2, 2, [1, 2]);
hold on;
plot(v, Laero);
plot(v, Daero);
legend("L", "D");
xlabel("Velocity (m/s)");
ylabel("Aerodynamic forces");
title("L and D function of velocity");
set(gca,'Xdir','reverse');
grid on;
axis([0 3700 0 190000]);

subplot(2, 2, 3)
hold on;
plot(v, Laero);
plot(v, Daero);
legend("L", "D");
xlabel("Velocity (m/s)");
ylabel("Aerodynamic forces");
title("L and D function of velocity while gliding");
set(gca,'Xdir','reverse');
grid on;
axis([vel_parachute 3700 0 80000]);

subplot(2, 2, 4);
hold on;
plot(v, Laero);
plot(v, Daero);
legend("L", "D");
xlabel("Velocity (m/s)");
ylabel("Aerodynamic forces");
title("L and D function of velocity the rest of the descent");
set(gca,'Xdir','reverse');
grid on;
axis([0 226.6 0 4000]);

%% gamma and gamma_ref function of time
figure;
% subplot(2, 2, [1, 2]);
grid on;
hold on;
plot(t, rad2deg(gamma));
plot(t, rad2deg(gamma_ref));
legend("\gamma", "\gamma_{ref}");
title("Flight path angle \gamma function of time");
xlabel("Time (s)");
ylabel("Flight path angle (deg)");
axis([0 t_landing -90 0]);
%%
% subplot(2, 2, 3);
grid on;
hold on;
plot(t, rad2deg(gamma));
plot(t, rad2deg(gamma_ref));
legend("\gamma", "\gamma_{ref}", 'Location','southwest');
title("Flight path angle \gamma function of time while gliding");
xlabel("Time (s)");
ylabel("Flight path angle (deg)");
axis([0 t_parachute -90 0]);

subplot(2, 2, 4);
grid on;
hold on;
plot(t, rad2deg(gamma));
plot(t, rad2deg(gamma_ref));
legend("\gamma", "\gamma_{ref}");
title("Flight path angle \gamma function of time the rest of the descent");
xlabel("Time (s)");
ylabel("Flight path angle (deg)");
axis([t_parachute t_landing -90 0]);

%% gamma and gamma_ref function of altitude
figure;
subplot(2, 2, [1, 2]);
grid on;
hold on;
plot(h, rad2deg(gamma));
plot(h, rad2deg(gamma_ref));
legend("\gamma", "\gamma_{ref}", 'Location','southwest');
title("Flight path angle \gamma function of altitude");
xlabel("Altitude (m)");
ylabel("Flight path angle (deg)");
axis([0 120000 -90 0]);
set(gca,'Xdir','reverse');

subplot(2, 2, 3);
grid on;
hold on;
plot(h, rad2deg(gamma));
plot(h, rad2deg(gamma_ref));
legend("\gamma", "\gamma_{ref}", 'Location','southwest');
title("Flight path angle \gamma function of altitude while gliding");
xlabel("Altitude (m)");
ylabel("Flight path angle (deg)");
axis([4000 120000 -90 0]);
set(gca,'Xdir','reverse');

subplot(2, 2, 4);
grid on;
hold on;
plot(h, rad2deg(gamma));
plot(h, rad2deg(gamma_ref));
legend("\gamma", "\gamma_{ref}");
title("Flight path angle \gamma function of altitude the rest of the descent");
xlabel("Altitude (m)");
ylabel("Flight path angle (deg)");
axis([0 4000 -90 0]);
set(gca,'Xdir','reverse');

%% Theta and alpha function of time
figure;
grid on;
hold on;
plot(t, rad2deg(theta));
plot(t, rad2deg(alpha));
legend("\theta", "\alpha");
title("\alpha and \theta function of time");
xlabel("Time (s)");
ylabel("Angles (deg)");
axis([0 t_landing -rad2deg(4) rad2deg(3)]);

%% M function of time
figure;
subplot(2, 2, [1, 2]);
plot(t, M);
title("Mach number function of time");
xlabel("Time (s)");
ylabel("Mach number");
grid on;
axis([0 t_landing 0 20]);

subplot(2, 2, 3);
plot(t, M);
title("Mach number function of time while gliding");
xlabel("Time (s)");
ylabel("Mach number");
axis([0 t_parachute 0.8 20]);
grid on;

subplot(2, 2, 4);
plot(t, M);
title("Mach number function of time the rest of the descent");
xlabel("Time (s)");
ylabel("Mach number");
axis([t_parachute t_landing 0 0.8]);
grid on;

%% Mach function of altitude
figure;
subplot(2, 3, [1, 3]);
plot(h, M);
xlabel("Altitude (m)");
ylabel("Mach number");
title("Mach number function of altitude");
set(gca,'Xdir','reverse');
axis([0 120000 0 20]);
grid on;

subplot(2, 3, 4);
plot(h, M);
xlabel("Altitude (m)");
ylabel("Mach number");
title("Mach number function of altitude while gliding");
axis([4000 120000 0.8 20]);
set(gca,'Xdir','reverse');
grid on;

subplot(2, 3, 5);
plot(h, M);
xlabel("Altitude (m)");
ylabel("Mach number");
title("Mach number function of altitude with parachute");
axis([150 4000 0 0.8]);
set(gca,'Xdir','reverse');
grid on;

subplot(2, 3, 6);
plot(h, M);
xlabel("Altitude (m)");
ylabel("Mach number");
title("Mach number function of altitude with retrorockets");
axis([0 150 0 0.05]);
set(gca,'Xdir','reverse');
grid on;


%% Capsule mass function of altitude
figure;
subplot(2, 2, [1, 2]);
plot(h, mass_capsule);
title("Capsule mass function of altitude");
ylabel("Mass");
xlabel("Altitude");
set(gca,'Xdir','reverse');
axis([0 120000 1000 1550]);
grid on;

subplot(2, 2, 3);
plot(h, mass_capsule);
title("Capsule mass function of altitude before retrorockets");
ylabel("Mass");
xlabel("Altitude");
set(gca,'Xdir','reverse');
axis([150 120000 1100 1550]);
grid on;

subplot(2, 2, 4);
plot(h, mass_capsule);
title("Capsule mass function of altitude with retrorockets");
ylabel("Mass");
xlabel("Altitude");
set(gca,'Xdir','reverse');
axis([0 150 1000 1100]);
grid on;

%% Heat flux
figure
subplot(2, 1, 1);
plot(t, heatflux);
title("Heatflux function of time")
xlabel("Time (s)");
ylabel("Heatflux (kW/m^2)")
grid on;
axis([0 t_landing 0 1800]);

subplot(2, 1, 2);
plot(t, integrated_heatflux);
title("Integrated heatflux function of time")
xlabel("Time (s)");
ylabel("Heatflux (kWh/m^2)")
axis([0 t_landing 0 11]);
grid on;

%% Thrust function of altitude
figure
plot(h(index_rockets:length(h)), thrust(index_rockets:length(h)),...
     h(index_rockets:length(h)), mass_capsule(index_rockets:length(h)).*g(index_rockets:length(h)));
legend("Thrust", "m\timesg");
title("Thrust function of altitude");
xlabel("Altitude (m)");
ylabel("Thrust (N)");
set(gca,'Xdir','reverse');
grid on;
print(gcf,'Thrust.png', '-dpng', '-r300')

%% Mass function of altitude
figure
plot(mass_capsule(1528: length(h)), h(1528: length(h)))
title("Mass variation during lander descent")
xlabel("Mass (kg)");
ylabel("Altitude (m)");
grid on;
axis([960 1080 -10 110]);
set(gca,'Xdir','reverse');
print(gcf,'MassVariation.png', '-dpng', '-r300')

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
