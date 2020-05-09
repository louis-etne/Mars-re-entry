function C_A = CASimulink(M, alpha)

    %% Observation : 
    % - from 0 to mach 1.5 : C_A depend only on mach number
    % - from mach 1.5 to ... : C_A depend on the angle of attack
    
    if M >=4
        C_A = -0.0004557 * alpha.^2 + 0.001214 * alpha + 1.59;
    elseif M >= 1.5 && M < 4
        C_A = -0.0005411 * alpha.^2 + 0.003171 * alpha + 1.599;
    else
        if abs(alpha) >= 17.5
            C_A = 0.4649 * M + 0.802;
        elseif abs(alpha) <= 17.5 && abs(alpha) >= 12.5
            C_A = 0.4229 * M + 0.848;
        elseif abs(alpha) <= 12.5 && abs(alpha) >= 8.5
            C_A = 0.4122 * M + 0.86;
        elseif abs(alpha) <= 8.5 && abs(alpha) >= 2.5
            C_A = 0.391 * M + 0.874;
        else
            C_A = 0.3714 * M + 0.88;
        end
    end
end
