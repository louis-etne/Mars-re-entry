clc
clear 
close all

%NAME: ROBERT SCHWEINSBERG
%MATR.NR.: 5002728


%EVERYTHING IN SI-UNITS
%**************************************************************************
%Input Parameter
G = 6.67408*10^-11;
m_e = 5.972*10^24; %Mass earth in SI
m_sc = 4100; %Mass space craft
r_e = 6378 * 10^3;%Radius  of Earth
A_wing = 70; %Wing Area
A_proj = A_wing; %Projected Area of SC
c_d = 1;%Drag coefficient
L_D = 0.1:0.7:2;
h0 = 35786 * 10^3;% Starting Altitude
m = 0.1; %Nose radius(Preliminary)
rn = 2.5; % 10*Radius of S/C(Preliminary)
r_proj = rn/2; % Radius of projected Area
nose_angle = 70;% Defines the Area of the nose sphere section
phi_max = nose_angle*((2*pi)/360); %Max Noseangle; %Max Noseangle
epsilon = 0.96; %Emission coefficient https://www.sika.net/images/Dokumente/Emissionsfaktoren_Tabelle.pdf
A_nose = 2*pi*rn^2*(1-cos(nose_angle*((2*pi)/360))); %Area of nose sphere section
v_sound = 295;%speed of sound


%**************************************************************************
%Atmospheric Parameters
BK = 9823.4;
XiK = 3.15;
PsiK = 0.5;
Br = 85.174;
Xir = 12.5;
Psir = 1.5;
rho0 = 1.226; %Atmospheric density at 0 m altitude

%**************************************************************************
%Initial Input Values
R0 = r_e + +h0; %Starting Radius
s0 = [R0 0 0]; %Input position vector
mue = m_e * G;
v_GEO = sqrt(mue/R0); %Velocity GEo
v_orbit = [0 v_GEO 0]; % Input velocoty vector
v0 = v_orbit -[0 1485.5 0]; %Deorbiting Vector

clear phi
syms phi

%**************************************************************************
%Heat shield
A_HS = A_wing; %Area of heatshield
m_spec_HSmaterial = 3210; %Siliziumcarb kg/m^3
d_HS = 0.003; %Thickness HS in m
V_HS = A_HS * d_HS; %Volume Heatshield
m_HS =  V_HS* m_spec_HSmaterial;%Mass heatshield
c_spec_HS = 750;  %Specific heat capacity J/(kg*K)


%**************************************************************************
%Solar Arrays 
E_solar = 1.367; %Solar konstant kW/m^2
W_max_burn = 1.6027*10^7;%Energy needed at GTO-GEO-burn
t_max_burn = 283; %Time of GTO-GEO burn
mue_solar_panel = 0.3; %Efficiency solar panels
t_max_cruise_recharge = 18560; %Cruise time GTO
Battery_efficiency = 0.85;% Battery efficiency
A_real = 0.8*A_wing; %realistic solar arry area
A_solar_panel = (W_max_burn)/(Battery_efficiency*mue_solar_panel*1000*E_solar*t_max_cruise_recharge); %Solar Panel area in m^2 for whole recharge for GTO_GEO burn
W_el = A_real*E_solar*1000*t_max_cruise_recharge*mue_solar_panel; %Electrical Energy that can be generated while cruise flight GTO 
rho_fuel_cell = 12/118;%Power density fuel cell in kW/kg ref: Space shuttle https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040010319.pdf
rho_SolarPanel = 0.8; %kW/kg 

%%
%Simulations & Output parameter
%**************************************************************************
c_a = 0.8;%Lift coefficient
sim('Atmospheric_Modell')
%Saving Output Data
s_1 = simout_s_.data; %Position vector array
alti1 = simout_h.data./1000; % Plotted in km
q1 = simout_q.data;%Accumulated heatflow
q_dot1 = simout_q_dot.data;%Heatflow
V1 = simout_V.data;%Velocity
inc1 = simout_i.data*360/(2* pi);%Inclination
a_ell1 = simout_a.data;%Semi major axis
e1 = simout_e.data;%Eccentricity
apogee1 = simout_apogee.data;%Apogee 
perigee1 = simout_perigee.data;%Perigee
h_am1 = simout_h_am.data;%Angular momentum
F_l1 = simout_F_l.data;%Lift
E_m1 = simout_spec_E_m.data;%Mechanical Energy (Orbit parameter)
T1 = simout_T.data;%temperatur heatshield(always facing sun)
acc1 = simout_acc.data;% G-forces 
t1 = simout_s_.time;%Time vector
SolarPanel_kW1 = simout_SolarPanel_kW.data;%Power output solar arrays
SolarPanel_kJ1 = simout_SolarPanel_kJ.data;%Energy output solar arrays
m_solarPanel1 = simout_SolarPanel_mass.data;%mass solar arrays
m_FuelCell1 = simout_FuelCell_mass.data;%Mass fuel cells without propellants
rho_1 = simout_rho.data;%Atmosphere density
T1_HS_shadow = simout_T_HS_shadow.data;%Temperature of heatshield always in shadow


%************************************************************************
c_a = 1.0;
sim('Atmospheric_Modell')
%Saving Output Data
s_2 = simout_s_.data;
alti2 = simout_h.data./1000; % Plotted in km
q2 = simout_q.data;
q_dot2 = simout_q_dot.data;
V2 = simout_V.data;
inc2 = simout_i.data*360/(2* pi);
a_ell2 = simout_a.data;
e2 = simout_e.data;
apogee2 = simout_apogee.data;
perigee2 = simout_perigee.data;
h_am2 = simout_h_am.data;
F_l2 = simout_F_l.data;
E_m2 = simout_spec_E_m.data;
T2 = simout_T.data;
acc2 = simout_acc.data;
t2 = simout_s_.time;
SolarPanel_kW2 = simout_SolarPanel_kW.data;
SolarPanel_kJ2 = simout_SolarPanel_kJ.data;
m_solarPanel2 = simout_SolarPanel_mass.data;
m_FuelCell2 = simout_FuelCell_mass.data;
T2_HS_shadow = simout_T_HS_shadow.data;

%************************************************************************
c_a = 1.5;
sim('Atmospheric_Modell')
%Saving Output Data
s_3 = simout_s_.data;
alti3 = simout_h.data./1000; % Plotted in km
q3 = simout_q.data;
q_dot3 = simout_q_dot.data;
V3 = simout_V.data;
inc3 = simout_i.data*360/(2* pi);
a_ell3 = simout_a.data;
e3 = simout_e.data;
apogee3 = simout_apogee.data;
perigee3 = simout_perigee.data;
h_am3 = simout_h_am.data;
F_l3 = simout_F_l.data;
E_m3 = simout_spec_E_m.data;
T3 = simout_T.data;
acc3 = simout_acc.data;
t3 = simout_s_.time;
SolarPanel_kW3 = simout_SolarPanel_kW.data;
SolarPanel_kJ3 = simout_SolarPanel_kJ.data;
m_solarPanel3 = simout_SolarPanel_mass.data;
m_FuelCell3 = simout_FuelCell_mass.data;
T3_HS_shadow = simout_T_HS_shadow.data;

%%
%Time in Atmosphere
n_points_in_atmosphere1 = find(alti1 < 200);
n_points_in_atmosphere2 = find(alti2 < 200);
n_points_in_atmosphere3 = find(alti3 < 200);

stepsize = 10; %Simulation intervalls in s

t_1_in_atm = length(n_points_in_atmosphere1)*stepsize; %Time in s
t_2_in_atm = length(n_points_in_atmosphere2)*stepsize; %Time in s
t_3_in_atm = length(n_points_in_atmosphere3)*stepsize; %Time in s

%%
%Heatshieldmass computation 
%PRELIMINARY COMPUTATION, NOT FURTHER EVALUATED
%ASSUMPTION: ALL HEAT INPUT IS ACTING ON THE NOSE AT 0° ANGLE OF ATTACK
%*************************************************************************
q_KWH1 = max(q1)*0.00027777777777778; % Accumulated Heat in KWH
q_KWH1 = q_KWH1 * cos(phi); % Lees-Law
    
m_spez_Ker1 = 0.32*q_KWH1 + 0.38*sqrt(q_KWH1) + 6.5;
m_spez_Metall1 = 0.34*q_KWH1+0.21*sqrt(q_KWH1) + 4.9;
m_spez_Ablator1 = 0.24*q_KWH1+0.29*sqrt(q_KWH1) + 11.3;
m_spez_water_active1 = 1.1 + 1.57 * q_KWH1;
m_spez_water_passive1 = 5.3 * 1.47* q_KWH1;

m_Ker1 = double(int(m_spez_Ker1 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_Metall1 = double(int(m_spez_Metall1 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_Ablator1 = double(int(m_spez_Ablator1 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_water_active1 = double(int(m_spez_water_active1 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_water_passive1 = double(int(m_spez_water_passive1 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));

%**************************************************************************
q_KWH2 = max(q2)*0.00027777777777778; % Accumulated Heat in KWH
q_KWH2 = q_KWH2 * cos(phi); % Lees-Law
    
m_spez_Ker2 = 0.32*q_KWH2 + 0.38*sqrt(q_KWH2) + 6.5;
m_spez_Metall2 = 0.34*q_KWH2+0.21*sqrt(q_KWH2) + 4.9;
m_spez_Ablator2 = 0.24*q_KWH2+0.29*sqrt(q_KWH2) + 11.3;
m_spez_water_active2 = 1.1 + 1.57 * q_KWH2;
m_spez_water_passive2 = 5.3 * 1.47* q_KWH2;

m_Ker2 = double(int(m_spez_Ker2 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_Metall2 = double(int(m_spez_Metall2 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_Ablator2 = double(int(m_spez_Ablator2 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_water_active2 = double(int(m_spez_water_active2 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_water_passive2 = double(int(m_spez_water_passive2 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));

%**************************************************************************
q_KWH3 = max(q3)*0.00027777777777778; % Accumulated Heat in KWH
q_KWH3 = q_KWH3 * cos(phi); % Lees-Law
    
m_spez_Ker3 = 0.32*q_KWH3 + 0.38*sqrt(q_KWH3) + 6.5;
m_spez_Metall3 = 0.34*q_KWH3+0.21*sqrt(q_KWH3) + 4.9;
m_spez_Ablator3 = 0.24*q_KWH3+0.29*sqrt(q_KWH3) + 11.3;
m_spez_water_active3 = 1.1 + 1.57 * q_KWH3;
m_spez_water_passive3 = 5.3 * 1.47* q_KWH3;

m_Ker3 = double(int(m_spez_Ker3 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_Metall3 = double(int(m_spez_Metall3 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_Ablator3 = double(int(m_spez_Ablator3 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_water_active3 = double(int(m_spez_water_active3 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));
m_water_passive3 = double(int(m_spez_water_passive3 * sin(phi)*rn^2 * 2*pi,phi,[0 phi_max]));

%**************************************************************************
%Putting Heat shield masses into arrays
m_metall = [m_Metall2; m_Metall2; m_Metall3];
m_keramics = [m_Ker1; m_Ker2; m_Ker3];
m_ablator = [m_Ablator1; m_Ablator2; m_Ablator3];
m_water_active = [m_water_active1; m_water_active2; m_water_active3];
m_water_passive = [m_water_passive1; m_water_passive2; m_water_passive3];

L_D_list = [0.8 1 1.5];

%%
%GRAPH PLOTTING

%Orbital Plot
figure('units','normalized','outerposition',[0 0 1 1])
original = imread('worldmap2.jpg');
[x y z] = ellipsoid(0,0,0,6378000,6378000,6378000);
surf(-x,-y,-z,'EdgeColor','none','LineStyle','none');
h_bla = findobj('Type','surface');
set(h_bla,'CData',original,'FaceColor','texturemap')
hold on
color = [0 0.5 0.3];
% plot3(s_1(:,1),s_1(:,2),s_1(:,3),'LineWidth',2,'Color',[1,0,0],'DisplayName','Orbit 1 [c_a = 0.8]');
% hold on 
% plot3(s_2(:,1),s_2(:,2),s_2(:,3),'LineWidth',2,'Color',[0,1,0],'DisplayName','Orbit 2 [c_a = 1.0]');
% hold on
plot3(s_3(:,1),s_3(:,2),s_3(:,3),'LineWidth',2,'Color',[0,0,1],'DisplayName','Orbit 3 [c_a = 1.5]');
legend
title('Atmospheric Breaking Manoeuvre Trajectory')
rectangle('Position',[-6778*10^3 -6778*10^3 6778*10^3*2 6778*10^3*2],'Curvature',[1 1],'Edgecolor',color, 'LineWidth',2) % LEO
rectangle('Position',[-42164*10^3 -42164*10^3 42164*10^3*2 42164*10^3*2],'Curvature',[1 1],'Edgecolor',[0.2 .0 .3], 'LineWidth',2) % GEO
rectangle('Position',[-6471*10^3 -6471*10^3 6471*10^3*2 6471*10^3*2],'Curvature',[1 1]) % Oberste Schicht
rectangle('Position',[-6382*10^3 -6382*10^3 6382*10^3*2 6382*10^3*2],'Curvature',[1 1]) % Grenze Troposphäre
rectangle('Position',[-6391*10^3 -6391*10^3 6391*10^3*2 6391*10^3*2],'Curvature',[1 1]) % Grenze Tropopause
rectangle('Position',[-6403*10^3 -6403*10^3 6403*10^3*2 6403*10^3*2],'Curvature',[1 1]) % Grenze Stratosphäre
ax = gca;
ax.Clipping = 'off';
axis equal

%**************************************************************************

%Other interesting Plots
% Altitude over time
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
% plot(t1,alti1,'r')
% hold on
% plot(t2,alti2,'g')
% hold on
plot(t3,alti3,'b')
xlabel('Time [s]')
ylabel('Altitude [km]')
grid on
legend('alti_1','alti_2','alti_3')
title('Altitude over time')

subplot(2,1,2)
% plot(t1,alti1,'r')
% hold on
% plot(t2,alti2,'g')
% hold on
plot(t3,alti3,'b')
xlabel('Time [s]')
ylabel('Altitude [km]')
ylim([105 125])
grid on
legend('alti_1','alti_2','alti_3')
title('Altitude over time')

% Velocity over Time
figure('units','normalized','outerposition',[0 0 1 1])
% plot(t1,V1,'r') 
% hold on
% plot(t2,V2,'g')
% hold on
plot(t3,V3,'b')
xlabel('Time [s]')
ylabel('Trajectory Velocity [m/s]')
grid on
legend('V_1','V_2','V_3')
title('Trajectory velocity over time')

% Velocity over Altitude
figure('units','normalized','outerposition',[0 0 1 1])
% plot(V1,alti1,'r') 
% hold on
% plot(V2,alti2,'g')
% hold on
plot(V3,alti3,'b')
xlabel('Trajectory Velocity [m/s]')
ylabel('Altitude [km]')
grid on
ylim([0 5000])
legend('V_1','V_2','V_3')
title('Trajectory velocity over altitude')

% Heatflow over time
figure('units','normalized','outerposition',[0 0 1 1])
% plot(t1,q_dot1,'r')
% hold on
% plot(t2,q_dot2,'g')
% hold on
plot(t3,q_dot3,'b')
xlabel('Time [s]')
ylabel('Heatflow [kW/m²]')
grid on
legend('q_d_o_t_1','q_d_o_t_2','q_d_o_t_3')
title('Heatflow over time')

% Accumulated Heatflow over Time
figure('units','normalized','outerposition',[0 0 1 1])
% plot(t1,q1,'r')
% hold on
% plot(t2,q2,'g')
% hold on
plot(t3,q3,'b')
xlabel('Time [s]')
ylabel('Accumulated Heatflow [kJ/m²]')
grid on
legend('q_1','q_2','q_3')
title('Accumulated heatflow over time')

%Eccentricity over time
figure('units','normalized','outerposition',[0 0 1 1])
% plot(t1,e1,'r')
% hold on
% plot(t2,e2,'g')
% hold on
plot(t3,e3,'b')
xlabel('Time [s]')
ylabel('eccentricity')
grid on
legend('e_1','e_2','e_3')
title('Eccentricity over time')

%Inclination over time
figure('units','normalized','outerposition',[0 0 1 1])
% plot(t1,inc1,'r')
% hold on
% plot(t2,inc2,'g')
% hold on
plot(t3,inc3,'b')
xlabel('Time [s]')
ylabel('Inclination [°]')
grid on
legend('inc_1','inc_2','inc_3')
title('Inclination over time')

%Apogee over time
figure('units','normalized','outerposition',[0 0 1 1])
plot(t1,apogee1,'r')
hold on
plot(t2,apogee2,'g')
hold on
plot(t3,apogee3,'b')
xlabel('Time [s]')
ylabel('Apogee [m]')
yline(400000+r_e);
grid on
legend('r_a_1','r_a_2','r_a_3')
title('Apogee over time')

%Spec. Mech. Energy over time
figure('units','normalized','outerposition',[0 0 1 1])
plot(t1,E_m1,'r')
hold on
plot(t2,E_m2,'g')
hold on
plot(t3,E_m3,'b')
xlabel('Time [s]')
ylabel('specific mech. energy [J]')
grid on
legend('E_m_1','E_m_2','E_m_3')
title('Spec. mech. energy over time')

%Semi-major-axis over time
figure('units','normalized','outerposition',[0 0 1 1])
plot(t1,a_ell1,'r')
hold on
plot(t2,a_ell2,'g')
hold on
plot(t3,a_ell3,'b')
xlabel('Time [s]')
ylabel('Semi-Major-axis [m]')
grid on
legend('a_e_l_l_1','a_e_l_l_2','a_e_l_l_3')
title('Semi-Major-Axis over time')

%Perigee over time
figure('units','normalized','outerposition',[0 0 1 1])
plot(t1,perigee1,'r')
hold on
plot(t2,perigee2,'g')
hold on
plot(t3,perigee3,'b')
xlabel('Time [s]')
ylabel('Perigee[m]')
grid on
legend('r_p_1','r_p_2','r_p_3')
title('Perigee over time')

%Plots Mass Heat shield
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(L_D_list,m_metall,'DisplayName','m_m_e_t_a_l')
hold on 
plot(L_D_list,m_keramics,'DisplayName','m_k_e_r_a_m_i_c_s')
hold on
plot(L_D_list,m_ablator,'DisplayName','m_a_b_l_a_t_o_r')
xlabel('L/D')
ylabel('Mass in [kg]')
legend
title('Heatshieldmass over L/D')

%Mass heatshields regression
subplot(2,1,2)
plot(L_D_list,m_water_active,'DisplayName','m_w_a_t_e_r_,_a_c_t_i_v_e')
hold on
plot(L_D_list,m_water_passive,'DisplayName','m_w_a_t_e_r_,_p_a_s_s_i_v_e')
xlabel('L/D')
ylabel('Mass in [kg]')
legend
title('Heatshieldmass(water cooled) over time')

%Temperature of heat shield
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
% plot(t1,T1,'r')
% hold on
% plot(t2,T2,'g')
% hold on
plot(t3,T3,'b')
yline(2300+273);
xlabel('t in s')
ylabel('T in K')
legend('T_1','T_2','T_3')
title('Temperature of heatshield (towards sun) over time')

subplot(2,1,2)
% plot(t1,T1_HS_shadow,'r')
% hold on
% plot(t2,T2_HS_shadow,'g')
% hold on
plot(t3,T3_HS_shadow,'b')
yline(2300+273);
xlabel('t in s')
ylabel('T in K')
legend('T_1_HS_shadow','T_2_HS_shadow','T_3_HS_shadow')
title('Temperature of heatshield (in shadow) over time')

%G_Forces
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
% plot(t1,acc1/9.81,'r')
% hold on
% plot(t2,acc2/9.81,'g')
% hold on
plot(t3,acc3/9.81,'b')
xlabel('t in s')
ylabel('g-forces')
legend('acceleration1','acceleration2','acceleration3')
title('acceleration over time')

%Acceleration
subplot(2,1,2)
% plot(t1,acc1,'r')
% hold on
% plot(t2,acc2,'g')
% hold on
plot(t3,acc3,'b')
xlabel('t in s')
ylabel('acceleration in m/s^2')
legend('acceleration1','acceleration2','acceleration3')
title('acceleration over time')

%Energy Solar Panels
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1)
plot(t1,SolarPanel_kW1,'r')
hold on
plot(t2,SolarPanel_kW2,'g')
hold on
plot(t3,SolarPanel_kW3,'b')
xlabel('t in s')
ylabel('Power in kW')
legend('E1','E2','E3')
title('Power over time')

subplot(2,1,2)
plot(t1,SolarPanel_kJ1,'r')
hold on
plot(t2,SolarPanel_kJ2,'g')
hold on
plot(t3,SolarPanel_kJ3,'b')
xlabel('t in s')
ylabel('Energy in kJ')
legend('E1','E2','E3')
title('Energy over time')

%Atmosphere
figure('units','normalized','outerposition',[0 0 1 1])
fplot(@(x_test) 1.226*exp(-x_test/8436),[0 210000])
xlabel('Altitude in m')
ylabel('Density in kg/m^3')
title('Atmospheric Model (Density over Altitude)')
