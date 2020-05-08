function f = Dynamics(t, y, Mars, Atm, Vehicle, aero_coefs, sim_params)
    % TODO
    % - CA and CN
    % - 
    
    % Assignations
    v = y(1);
    gamma = y(2);
    h = y(3);
    phi = y(4);
    theta = y(5);
    q = y(6);
    % m = y(5) To be implemented
    
    % Intermediate calculation
    r = h + Mars.radius * 1000;
    g = Mars.mu / r^2; % m/s^2 - Gravitationnal acceleration at r
    rho = Density(h, Atm); % kg/m^3 - Air density at h
    Pdyn = (1/2) * rho * v^2; % Dynamic pressure
    M = MachNumber(v, h, Atm); % Mach number
    alpha = theta - gamma;
    
    CA = AxialForceCoef(M, alpha, aero_coefs.CA);
    CN = NormalForceCoef(M, alpha, aero_coefs.CN);
    
    CL = CN;
    CD = CA;
    
    if (sim_params.controlled)
        B = CD * Vehicle.S / Vehicle.mass;
        dv_aero = sim_params.v_des - sqrt(v^2+((2*Mars.mu)*((1/sim_params.r_des)-(1/r))));
        gamma_ref = asin((1/2)*B*Atm.hs*((sim_params.rho_des - rho) ./ log(1 + dv_aero/v)));
        
        k = (Pdyn * Vehicle.S * CL) / (v * Vehicle.mass); % equivalent to g
        theta_eq = gamma - (cos(gamma) * Vehicle.mass) / (Pdyn*Vehicle.S*CL) * (v^2/r - Mars.mu / r^2);
        theta_cmd = theta_eq + (1/(sim_params.tau*k)) * (gamma_ref - gamma);
        
        % saturation
        if (theta_cmd > deg2rad(60))
            theta_cmd = deg2rad(60);
        end
        
        Kp = sim_params.wn^2;
        Kd = 2 * sim_params.zeta * sim_params.wn;
        k = (Vehicle.CMdelta * Pdyn * Vehicle.S * Vehicle.d) / Vehicle.J;
        delta_eq = (-Vehicle.CMalpha/Vehicle.CMdelta) * (theta - gamma) - (Vehicle.CMq * Vehicle.d * q) / (2 * v * Vehicle.CMdelta);
        delta_cmd = delta_eq + Kp/k * (theta_cmd - theta) + Kd/k * (-q);
    else
        delta_cmd = 0;
    end
    
    Daero = Pdyn * Vehicle.S * CD;
    Laero = Pdyn * Vehicle.S * CL * alpha;
    Maero = Pdyn * Vehicle.S * Vehicle.d * (Vehicle.CMalpha * alpha + (Vehicle.d / (2*v)) * Vehicle.CMq*q + Vehicle.CMdelta * delta_cmd);
    
    if (h <= sim_params.parachute_at)
        Daero = Daero + Vehicle.Parachute.CD .* ((rho .* v.^2) / 2) * Vehicle.Parachute.S;
    end
    
    if (h <= 0)
        f(1) = 0;
    else
        f(1) = -Daero / Vehicle.mass - g*sin(gamma);
    end
    
    f(2) = (1/v) * (Laero/Vehicle.mass + ((v.^2/r) - g)*cos(gamma));
    
    if (h <= 0)
        f(3) = 0;
    else
        f(3) = v * sin(gamma);
    end
    f(4) = v/r * cos(gamma);
    f(5) = q;
    f(6) = (1 / Vehicle.J) * Maero;
    
    f = f(:);
end

