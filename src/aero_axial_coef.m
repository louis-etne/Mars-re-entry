
alpha = 0:25;
C_A_M6_10 = [1.6 1.6 1.59 1.59 1.59 1.58 1.57 1.57 1.565 1.56 1.55 1.54 1.535 1.53 ...
       1.52 1.51 1.50 1.495 1.48 1.45 1.43 1.42 1.40 1.38 1.35 1.32];
   
figure();
plot(alpha, C_A);

f_CA_6_10 = fit(transpose(alpha), transpose(C_A), 'poly2');

alpha_test = 0:0.01:25
c_test = f_CA.p1 * alpha_test.^2 + f_CA.p2 * alpha_test + f_CA.p3

hold on;
plot(alpha_test, c_test);
