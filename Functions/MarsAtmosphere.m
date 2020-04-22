function [T, P, rho] = MarsAtmosphere(h)
    function t = Temperature(hi)
        if hi > 7000
            t = -23.4 - 0.00222 * hi; % °C - Temperature
        else
            t = -31.0 - 0.000998 * hi; % °C - Temperature
        end
    end

    function rho = Density(h)
        addpath('Data');
        load('constants.mat', 'Atm');

        rho = Atm.rho0 * exp(-h/Atm.hs);
    end

T = arrayfun(@Temperature, h);
P = 0.699 .* exp(-0.00009 .* h);

% Calculation of density
rho = Density(h); % kg/m^3
end