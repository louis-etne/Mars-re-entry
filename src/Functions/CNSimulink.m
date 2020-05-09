function C_N = CNSimulink(M, alpha)
    if M >= 8
        C_N = 0.00282 * alpha + -0.003312;
    elseif M >= 4 && M <= 8 
        C_N = 0.002487 * alpha + -0.0009416;
    elseif M >= 1.5 && M <= 4
        C_N = -9.821e-05 * alpha.^2 + -0.002146 * alpha + -0.0002292;
    elseif M >= 1.1 && M <= 1.5
        C_N = 2.432e-06 * alpha.^2 + 0.00231 * alpha + 0.01058;
    elseif M >= 0.95 && M <= 1.1 
        C_N = -3.433e-05 * alpha.^2 + 0.001611 * alpha + 0.007388;
    elseif M >= 0.7 && M <= 0.95 
        C_N = -4.826e-05 * alpha.^2 + 0.001155 * alpha + 0.004005;
    elseif M >= 0.5 && M <= 0.7 
        C_N = -7.615e-05 * alpha.^2 + 0.000277 * alpha + 0.00178;
    else 
        if abs(alpha) >=5
            C_N = -4.718e-05 * alpha.^2 + 0.0008791 * alpha + 0.0056;
        elseif abs(alpha) <= 1
            C_N = 0.002 * alpha + 0.003;
        else
            C_N = 0;
        end
    end
end