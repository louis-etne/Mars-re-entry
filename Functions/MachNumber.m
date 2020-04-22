function [mach] = MachNumber(T)
    addpath('Data');
    load('constants.mat', 'Atm');

    mach = sqrt(Atm.gamma * Atm.R * T);
end