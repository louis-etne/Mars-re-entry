function [T, P, rho] = mars_atmos(h)
if h > 7000 %check alt
    T = -23.4-0.00222 * h; %Temp in celcius
    P = 0.699 * exp(-0.00009 * h); % Pressure (K -Pa)
else
    T = -31 - 0.000998 * h;
    P = 0.699 * exp(-0.00009 * h);
end

% Calculation of density
rho = P / (0.1921 * (T + 273.1)); %kg/m^3
end

