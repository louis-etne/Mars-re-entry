clear all;
close all;
clc;

% NASA model stops at 100km
h1 = 0:2:100;
% Teacher models stop at 150km but is less precise
h2 = 110:10:150;
% So we use the NASA model first and after the teacher
h = [h1 h2] * 1000; % Convert h in meters

% Temperature is in K
t1 = [214, 213.8, 213.4, 212.4, 209.2, 205, 201.4, 197.8,...
      194.6, 191.4, 188.2, 185.2, 182.5, 180, 177.5, 175, 172.5, 170,...
      167.5, 164.8, 162.4, 160, 158, 156, 154.1, 152.2, 150.3, 148.7,...
      147.2, 145.7, 144.2, 143, 142, 141, 140, 139.5, 139, 139, 139,...
      139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139];
t2 = [149.4, 159.7, 170, 245.1, 288.6];
T = [t1 t2];

plot(T, h);

%% Fits models model

%from 0 to 10 : 
T_0_10 = t1(1 : 6);
h_0_10 = (0:2:10) * 1e3;

f_0_10 = fit(transpose(h_0_10), transpose(T_0_10), 'poly3');

h_test = (0:0.1:10) * 1e3;
t_test = f_0_10.p1 * h_test.^3 + f_0_10.p2 * h_test.^2 + ...
         f_0_10.p3 * h_test.^1 + f_0_10.p4;

figure();
plot(h_0_10, T_0_10);
hold on;
plot(h_test, t_test);

% from 10 to 70 : 
T_10_70 = t1(6 : 36);
h_10_70 = (10:2:70) * 1e3;
f_10_70=fit(transpose(h_10_70), transpose(T_10_70), 'poly4');

h_test = (10:0.1:70) * 1e3;
t_test = f_10_70.p1 * h_test.^4 + f_10_70.p2 * h_test.^3 + ...
     f_10_70.p3 * h_test.^2 + f_10_70.p4 * h_test.^1 + f_10_70.p5;

figure();
plot(h_10_70, T_10_70);
hold on;
plot(h_test, t_test);
 
% from 70 to 100 : line at 139

%from 100 to 150
 
h2 = (100 : 10 : 120) * 1e3;
t2 = [139, 149.4, 159.7];
f_100_150 = fit(transpose(h2), transpose(t2), 'poly1');

h_test_100_150 = (100:0.1:120) * 1e3;
t_test_100_150 = f_100_150.p1 * h_test_100_150 + f_100_150.p2;

figure();
plot(h2, t2);
hold on;
plot(h_test_100_150, t_test_100_150);

save('Data\fit_temp', 'f_0_10', 'f_10_70', 'f_100_150');

% compare all;
clear f_0_10 f_10_70 f_100_150;
addpath('Functions');
figure();
plot(T, h);
hold on;
plot(Temperature(h), h);
