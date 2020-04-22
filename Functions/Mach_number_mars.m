function [mach] = Mach_number_mars(T)
gamma_mars = 1.29;
R_mars = 191.8;
mach = sqrt(gamma_mars * R_mars * T)
end

