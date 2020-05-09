function C_A = AxialForceCoef(M, alpha, coefs)

    %% Observation : 
    % - from 0 to mach 1.5 : C_A depend only on mach number
    % - from mach 1.5 to ... : C_A depend on the angle of attack
    
    if M >=4
        C_A = coefs.f_CA_6_10.p1 * alpha.^2 + coefs.f_CA_6_10.p2 * alpha + coefs.f_CA_6_10.p3;
    elseif M >= 1.5 && M < 4
        C_A = coefs.f_CA_2.p1 * alpha.^2 + coefs.f_CA_2.p2 * alpha + coefs.f_CA_2.p3;
    else
        if abs(alpha) >= 17.5
            C_A = coefs.f_CA_a20.p1 * M + coefs.f_CA_a20.p2;
        elseif abs(alpha) <= 17.5 && abs(alpha) >= 12.5
            C_A = coefs.f_CA_a15.p1 * M + coefs.f_CA_a15.p2;
        elseif abs(alpha) <= 12.5 && abs(alpha) >= 8.5
            C_A = coefs.f_CA_a10.p1 * M + coefs.f_CA_a10.p2;
        elseif abs(alpha) <= 8.5 && abs(alpha) >= 2.5
            C_A = coefs.f_CA_a5.p1 * M + coefs.f_CA_a5.p2;
        else
            C_A = coefs.f_CA_a0.p1 * M + coefs.f_CA_a0.p2;
        end
    end
end