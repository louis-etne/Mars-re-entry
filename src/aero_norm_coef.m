%% Mach 6
alpha = 0 : 20;

C_L_M6  = [0 0.002 0.005 0.006 0.009 0.011 0.014 0.016 0.019 0.021 ...
           0.023 0.026 0.028 0.032 0.034 0.036 0.0385 0.041 0.044 0.047 0.05];
      
f_CL6 = fit(transpose(alpha), transpose(C_L_M6),  'poly1');

alpha_test = 0 : 0.01 : 20;
C_test = f_CL6.p1 * alpha_test + f_CL6.p2;

figure();
plot(alpha, C_L_M6);
hold on;
plot(alpha_test, C_test);

%% Mach 10
alpha = 0 : 25;

C_L_M10 = [0 0.002 0.005 0.006 0.009 0.011 0.014 0.016 0.019 0.021 ...
           0.023 0.026 0.028 0.032 0.034 0.036 0.0385 0.042 0.045 0.048 0.053 ...
           0.056 0.06 0.063 0.068 0.075];
       
f_CL10 = fit(transpose(alpha), transpose(C_L_M10), 'poly1');

alpha_test = 0 : 0.01 : 20;
C_test = f_CL10.p1 * alpha_test + f_CL10.p2;
       
figure();
plot(alpha, C_L_M10);
hold on;
plot(alpha_test, C_test);

    
addpath('Data');
save('Data\fit_areo_norm_coef', 'C_L_M6', 'C_L_M10');