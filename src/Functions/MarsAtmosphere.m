function [T, P, rho] = MarsAtmosphere(h)
    T = Temperature(h); % K - Temperature function of h
    P = Pressure(h); % Pascal - Pressure function of h
    rho = Density(h); % kg/m^3 - Density function of h
end