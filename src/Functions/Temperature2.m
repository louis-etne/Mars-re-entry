function T = Temperature2(h)    
    
    load('Data\fit_temp', 'f_0_10', 'f_10_70', 'f_100_150');

    function t = temp(hi)
        if hi < 10000
            t = f_0_10.p1 * hi.^3 + f_0_10.p2 * hi.^2 +...
                f_0_10.p3 * hi.^1 + f_0_10.p4; % K - Temperature
        elseif hi >= 10000 && hi <= 70000
            t = f_10_70.p1 * hi.^4 + f_10_70.p2 * hi.^3 + ...
                  f_10_70.p3 * hi.^2 + f_10_70.p4 * hi.^1 + ...
                  f_10_70.p5; % K - Temperature
        elseif hi > 70000 && hi <= 100000
            t = 139; % K - Temperature
        else
            t = f_100_150.a * exp(hi * f_100_150.b) +...
                f_100_150.c * exp(hi * f_100_150.d); % K - Temperature
        end
    end

    T = arrayfun(@temp, h); % K - Temperature
end
