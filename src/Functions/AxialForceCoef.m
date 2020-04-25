function C_A = AxialForceCoef(M, alpha)

    %% Observation : 
    % - from 0 to mach 1.5 : C_A depend only on mach number
    % - from mach 1.5 to ... : C_A depend on the angle of attack
    addpath('Data');
    load('fit_areo_axial_coef.mat', 'f_CA_2', 'f_CA_6_10', ...
        'f_CA_a20', 'f_CA_a15', 'f_CA_a10', 'f_CA_a5', 'f_CA_a0');

    if M >=4
        C_A = f_CA_6_10.p1 * alpha.^2 + f_CA_6_10.p2 * alpha + f_CA_6_10.p3;
    elseif M >= 1.5 && M < 4
        C_A = f_CA_2.p1 * alpha.^2 + f_CA_2.p2 * alpha + f_CA_2.p3;
    elseif M < 1.5
        if abs(alpha) >= 17.5
            C_A = f_CA_a20.p1 * M + f_CA_a20.p2;
        elseif abs(alpha) <= 17.5 && abs(alpha) >= 12.5
            C_A = f_CA_a15.p1 * M + f_CA_a15.p2;
        elseif abs(alpha) <= 12.5 && abs(alpha) >= 8.5
            C_A = f_CA_a15.p1 * M + f_CA_a15.p2;
        elseif abs(alpha) <= 8.5 && abs(alpha) >= 2.5
            C_A = f_CA_a15.p1 * M + f_CA_a15.p2;
        else
            C_A = f_CA_a15.p1 * M + f_CA_a15.p2;
        end
    end
    
end