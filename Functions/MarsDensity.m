function rho = MarsDensity(h)
    addpath('Data');
    load('constants.mat', 'hs', 'rho0');
    
    rho = rho0 * exp(-h/hs);
end

