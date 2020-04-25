function C_N = NormalForceCoef(normal_acc, vel, h)
    addpath('Data');
    load('constants.mat', 'Vehicle');
    
    rho = Density(h);
    
    C_N = (Vehicle.mass .* normal_acc) ./ (0.5 .* rho .* vel.^2 .* Vehicle.S);
end