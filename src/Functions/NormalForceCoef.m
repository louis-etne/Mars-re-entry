function C_N = NormalForceCoef(M, alpha, coefs)
    if M >= 8
        C_N = coefs.f_CN10.p1 * alpha + coefs.f_CN10.p2;
    elseif M >= 4 && M <= 8 
        C_N = coefs.f_CN10.p1 * alpha + coefs.f_CN10.p2;
    elseif M >= 1.5 && M <= 4
        C_N = coefs.f_CN_2.p1 * alpha.^2 + coefs.f_CN_2.p2 * alpha + coefs.f_CN_2.p3;
    elseif M >= 1.1 && M <= 1.5
        C_N = coefs.f_CN_1_2.p1 * alpha.^2 + coefs.f_CN_1_2.p2 * alpha + coefs.f_CN_1_2.p3;
    elseif M >= 0.95 && M <= 1.1 
        C_N = coefs.f_CN_1.p1 * alpha.^2 + coefs.f_CN_1.p2 * alpha + coefs.f_CN_1.p3;
    elseif M >= 0.7 && M <= 0.95 
        C_N = coefs.f_CN_08.p1 * alpha.^2 + coefs.f_CN_08.p2 * alpha + coefs.f_CN_08.p3;
    elseif M >= 0.5 && M <= 0.7 
        C_N = coefs.f_CN_06.p1 * alpha.^2 + coefs.f_CN_06.p2 * alpha + coefs.f_CN_06.p3;
    else 
        if abs(alpha) >=5
            C_N = coefs.f_CN_04_a_23_5.p1 * alpha.^2 + coefs.f_CN_04_a_23_5.p2 * alpha ...
                + coefs.f_CN_04_a_23_5.p3;
        elseif abs(alpha) <= 1
            C_N = coefs.f_CN_04_a_1_0.p1 * alpha + coefs.f_CN_04_a_1_0.p2;
        else
            C_N = 0;
        end
    end
end