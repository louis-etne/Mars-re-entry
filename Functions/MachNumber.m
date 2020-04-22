function [mach] = MachNumber(T)
    addpath('Data');
    load('constants.mat', 'gamma_m', 'R_m');

    mach = sqrt(gamma_m * R_m * T);
end