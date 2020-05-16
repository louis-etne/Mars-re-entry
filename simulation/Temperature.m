function T = Temperature(h)
    function t = temp(hi)
        if hi < 10000
            t = -1.019e-11 * hi.^3 + 1.349e-08 * hi.^2 +...
                -1.958e-05 * hi.^1 + 214; % K - Temperature
        elseif hi >= 10000 && hi <= 70000
            t = 2.845e-18 * hi.^4 + -4.139e-13 * hi.^3 + ...
                  3.048e-08 * hi.^2 -0.002304 * hi.^1 + ...
                  225.2; % K - Temperature
        elseif hi > 70000 && hi <= 100000
            t = 139; % K - Temperature
        else
            t = 0.001035 * hi + 35.52; % K - Temperature
        end
    end

    T = arrayfun(@temp, h); % K - Temperature
end
