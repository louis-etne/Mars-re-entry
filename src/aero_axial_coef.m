close all;
clear all;
clc;

alpha = 0:25;
C_A_M6_10 = [1.6 1.6 1.59 1.59 1.59 1.58 1.57 1.57 1.565 1.56 1.55 1.54 1.535 1.53 ...
       1.52 1.51 1.50 1.495 1.48 1.45 1.43 1.42 1.40 1.38 1.35 1.32];
   
figure();
plot(alpha, C_A_M6_10);

f_CA_6_10 = fit(transpose(alpha), transpose(C_A_M6_10), 'poly2');

alpha_test = -25:0.01:25;
c_test = f_CA_6_10.p1 * alpha_test.^2 + f_CA_6_10.p2 * alpha_test + f_CA_6_10.p3;

hold on;
plot(alpha_test, c_test);

%% M = 2
alpha = -8:0;
C_A_M2 = [1.54 1.55 1.56 1.57 1.58 1.585 1.59 1.595 1.6 ];

f_CA_2 = fit(transpose(alpha), transpose(C_A_M2), 'poly2');

alpha_test = -25:0.01:25;
c_test2 = f_CA_2.p1 * alpha_test.^2 + f_CA_2.p2 * alpha_test + f_CA_2.p3;

figure();
plot(alpha, C_A_M2);
ylim([1 1.8])
xlim([-25 25])
hold on;
plot(alpha_test, c_test2);

figure();
plot(alpha_test, c_test);
hold on;
plot(alpha_test, c_test2);

%% for other M 

M = [0.4 0.6 0.8 0.9 1.00 1.20];
a20 = [1.01 1.07 1.13 1.2 1.33 1.35];
a15 = [1.03 1.09 1.17 1.2 1.33 1.34];
a10 = [1.03 1.09 1.19 1.22 1.32 1.33];
a5 = [1.04 1.09 1.18 1.21 1.32 1.32];
a0 = [1.04 1.09 1.16 1.2 1.3 1.31];

f_CA_a20 = fit(transpose(M), transpose(a20), 'poly1');
f_CA_a15 = fit(transpose(M), transpose(a20), 'poly1');
f_CA_a10 = fit(transpose(M), transpose(a20), 'poly1');
f_CA_a5 = fit(transpose(M), transpose(a20), 'poly1');
f_CA_a0 = fit(transpose(M), transpose(a20), 'poly1');

addpath('Data');
save('Data/fit_areo_axial_coef.mat', 'f_CA_6_10', 'f_CA_2', ...
    'f_CA_a20', 'f_CA_a15', 'f_CA_a10', 'f_CA_a5', 'f_CA_a0');
