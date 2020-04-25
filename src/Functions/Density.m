function rho = Density(h)
    addpath('Data');
    load('constants.mat', 'Atm');
    
    function r = tmp(hi)
        if hi >= 120000
            r = 0;
        else
            r =  Atm.rho0 * exp(-hi/Atm.hs);
        end
    end
    
    rho = arrayfun(@tmp, h);
end