function rho = Density(h)
    addpath('Data');
    load('constants.mat', 'Atm');

    rho = Atm.rho0 * exp(-h/Atm.hs);
end