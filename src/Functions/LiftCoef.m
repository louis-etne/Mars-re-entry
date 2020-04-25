function C_L = LiftCoef(h, alpha, vel, axial_acc, normal_acc)
    % h : m - altitude
    % alpha: rad - angle of attack
    % vel: m/s - velocity
    % axial_acc : m/s^2 - axial acceleration
    % normal_acc : m/s^2 - normal acceleration
    
    C_A = AxialForceCoef(axial_acc, vel, h);
    C_N = NormalForceCoef(normal_acc, vel, h);
    
    C_L = -C_A * sin(alpha) + C_N * cos(alpha);
end