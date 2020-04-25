function C_A = AxialForceCoef(axial_acc, vel, h)
    addpath('Data');
    load('constants.mat', 'Vehicle');
    
    rho = Density(h);
    
    C_A = (Vehicle.mass .* axial_acc) ./ (0.5 .* rho .* vel.^2 .* Vehicle.S);
end