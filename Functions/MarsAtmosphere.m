function [T, P, rho1, rho2] = MarsAtmosphere(h)
    function t = Temperature(hi)
        if hi > 7000
            t = -23.4 - 0.00222 * hi; % °C - Temperature
        else
            t = -31.0 - 0.000998 * hi; % °C - Temperature
        end
    end

    function rho = Density(h)
        addpath('Data');
        load('constants.mat', 'hs', 'rho0');

        rho = rho0 * exp(-h/hs);
    end

T = arrayfun(@Temperature, h);
P = 0.699 .* exp(-0.00009 .* h);

% Calculation of density
rho1 = Density(h); %kg/m^3
rho2 = P ./ (0.1921 .* (T + 273.15));
end