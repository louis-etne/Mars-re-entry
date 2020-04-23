function T = Temperature(h)
    function t = temp(hi)
        if hi > 7000
            t = -23.4 - 0.00222 * hi; % °C - Temperature
        else
            t = -31.0 - 0.000998 * hi; % °C - Temperature
        end
    end

    T = arrayfun(@temp, h) + 273.15; % K - Temperature
end

