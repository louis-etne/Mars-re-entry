function P = Pressure(h)
    function p = tmp(hi)
        if hi >= 1200000
            p = 0;
        else
            p = 0.699 .* exp(-0.00009 .* hi);
        end
    end
    
    P = arrayfun(@tmp, h);
end

