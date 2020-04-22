%% Parameter of orbit around mars : 
% HYPOTHESIS : 
%   - circular orbit at 500 km altitude

%% constant :

alt = 500; %km
radius_m = 3389.5;
mu_m = 42828;
% circular orbit : 
r = radius_m + alt;

% speed departure orbit
V_departure_orbit = sqrt(2 * ( -mu_m / (2 * r) + mu_m / r));

%% 