function f = Dynamics(t, y, Mars, Atm, Vehicle)
    % Don't import data with load or use function in Dynamics,
    % It's really slow otherwise

    % Assignations
    h = y(1);
    v = y(2);
    s = y(3);
    gamma = y(4);
    theta = y(5);
    q = y(6);
    
    % Intermediate calculation
    r = Mars.radius * 1000 + h;
    g = Mars.mu / r^2; % m/s^2 - Gravitationnal acceleration at r
    alpha = theta - gamma; % rad - Angle of attack
    rho = Atm.rho0 * exp(-h/Atm.hs); % kg/m^3 - Air density at h
    Pdyn = (1/2) * rho * v^2; % Dynamic pressure
    
    Daero = Pdyn * Vehicle.S * Vehicle.CD0;
    Laero = Pdyn * Vehicle.S * Vehicle.CLalpha * alpha;
    Maero = Pdyn * Vehicle.S * Vehicle.d * (Vehicle.CMalpha * alpha + (Vehicle.d / (2*v)) * Vehicle.CMq*q);
    
    if (h <= 0)
        f(1) = 0;
    else
        f(1) = v*sin(gamma);
    end
    f(2) = -Daero / Vehicle.mass - g*sin(gamma);
    f(3) = (v/r) * cos(gamma);
    f(4) = (1/v) * (Laero/Vehicle.mass + ((v.^2/r) - g)*cos(gamma));
    f(5) = q;
    f(6) = (1/Vehicle.J) * Maero;
    
    f = f(:);
end

