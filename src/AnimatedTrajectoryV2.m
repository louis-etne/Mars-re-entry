clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat', 'Mars', 'Orbit', 'Vehicle');
load('altitude.mat', 'alt', 'downrange', 'time');
%% Parameter of orbit around mars : 
% HYPOTHESIS : 
%   - circular orbit at 500 km altitude
% Supposition for the calculation : Mars Atmosphere begin at 120 km high
% speed departure orbit
V_departure_orbit = sqrt(2 * (-Mars.mu / (2 * Orbit.radius) + Mars.mu / Orbit.radius));

r_p = Mars.radius + 0;   % km - altitude of transfer orbit periapsis
a_i = Orbit.radius;      % km - semi-major axis of initial orbit
e_i = 0;                 % eccentricity of initial orbit

a_t = (Orbit.radius + r_p)/2;                        % km - sma of transfer orbit
e_t = (Orbit.radius - r_p) / (Orbit.radius + r_p);   % eccentricity of transfer orbit
Va_t = sqrt(Mars.mu * (2/Orbit.radius - 1/a_t));     % km/s velocity at apoapsis of transfer orbit
Vp_t = sqrt(Mars.mu * (2/r_p - 1/a_t));              % km/s velocity at periapsis of transfer orbit
DeltaV = Va_t - V_departure_orbit;
h = r_p * Vp_t;
p = h^2 / Mars.mu;
r_entry = Mars.radius + 120;
theta_entry = acos(p/(r_entry * e_t) - 1/e_t);
V_entry = sqrt(Mars.mu * (2/r_entry - 1/a_t));

InitialOrbit = flip(Orbites(a_i*1e3, e_i, 0, 0)./1e3);     % divide by 1e3 because orbites fucntion returns values in meters
TransferOrbit = flip(Orbites(a_t*1e3, e_t, 0, 0)./1e3);

lon = rad2deg(theta_entry) - rad2deg(downrange);
lat = zeros(length(alt), 1);
lat = lat + 1;
alt = (alt./1e3 + Mars.radius);
Pos = SphereToCartesian(alt, lon, 90);
ellipse1 = animatedline('Color', 'b', 'LineWidth', 1.5);
ellipse2 = animatedline('Color', 'c', 'LineWidth', 1.5);
descent = animatedline('Color', 'g', 'LineWidth', 2);

%% Input Handling
% [cax,args,nargs] = axescheck(varargin{:}); % Parse possible Axes input
% error(nargchk(0,2,nargs)); % Ensure there are a valid number of inputs
% Handle remaining inputs.
% Should have 0 or 1 string input, 0 or 1 numeric input
j = 0;
k = 0;
n = 50; % default value
units = 'km'; % default value
% for i = 1:nargs
%     if ischar(args{i})
%         units = args{i};
%         j = j+1;
%     elseif isnumeric(args{i})
%         n = args{i};
%         k = k+1;
%     end
% end
if j > 1 || k > 1
    error('Invalid input types')
end
%% Calculations
% Scale factors
Scale = {'km' 'm'  'mile'            'miles'           'nm'              'au'                 'ft';
         1    1000 0.621371192237334 0.621371192237334 0.539956803455724 6.6845871226706e-009 3280.839895};
% Identify which scale to use
try
    myscale = Mars.radius*Scale{2,strcmpi(Scale(1,:),units)};
catch %#ok<*CTCH>
    error('Invalid units requested. Please use m, km, ft, mile, miles, nm, or AU')
end
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = myscale*cosphi*cos(theta);
y = myscale*cosphi*sintheta;
z = myscale*sin(phi)*ones(1,n+1);
%% Plotting
%     cax = newplot(cax);
    % Load and define topographic data
%     figure('units','normalized','outerposition',[0 0 1 1])
    load('topo.mat','topo','topomap1');
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    hold on
    topo2 = [topo(:,181:360) topo(:,1:180)]; %# ok<NODEF>
    mars_texture = imread('./Data/2k_mars.jpg');
    % Define surface settings
    props.Cdata = mars_texture;
    props.FaceColor= 'texturemap';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    % Create the sphere with Earth topography and adjust colormap
%     surface(x,y,z,props,'parent',gca)
    colormap(topomap1)
    
    axlen = Mars.radius*1.2;
    axis([-axlen, axlen, -axlen, axlen, -axlen, axlen])
%     axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])

    xx = x; yy = y; zz = z;

    ax = gca;
    ax.Clipping = 'off';
    
    view(127.5,30)
    pause(1);
    index1 = int16(length(InitialOrbit)/2);
    index2 = int16(length(TransferOrbit)/2);
    i1 = 1;
    i2 = index2;
    i3 = 1;
    for t = linspace(0,0.1,1000)
            % ----------    PLAY WITH THIS
%         st1 = sind(sind(t) );
%         ct1 = cosd(sind(t) );
        st3 = sind(t);
        ct3 = cosd(t);
            % ------------
        Rx = [1 0 0; 0 1 0; 0 0 1];    % rotation matrix around X axis
        Rz = [ct3 -st3 0; st3 ct3 0; 0 0 1];    % rotation matrix around Z axis
        V = Rx*[x(:) y(:) z(:)]';               % rotate data
        V = Rz*V;
        x = reshape(V(1,:),size(x));
        y = reshape(V(2,:),size(y));
        z = reshape(V(3,:),size(z));
        h = surface(x, y, z, props,'parent',gca);
        pause(0.005)
        if i1 <= index1
            addpoints(ellipse1, InitialOrbit(i1, 1), InitialOrbit(i1, 2), InitialOrbit(i1, 3));
            drawnow;
            i1 = i1+1;
        end
        if i1 > index1 && i2 <= 302
            addpoints(ellipse2, TransferOrbit(i2, 1), TransferOrbit(i2, 2), TransferOrbit(i2, 3));
            drawnow;
            i2 = i2+1;
        end
        if i2 > 302 && i3 <= length(Pos)
            addpoints(descent, Pos(i3, 1), Pos(i3, 2), Pos(i3, 3));
            drawnow;
            i3 = i3 + 1;
        end
        if t ~= 0.1
            delete(h)
        end
    end


%%

hold on


% for i = 1:int16(length(InitialOrbit)/2)
%    addpoints(ellipse1, InitialOrbit(i, 1), InitialOrbit(i, 2), InitialOrbit(i, 3));
%    drawnow;
% end
% for i = int16(length(TransferOrbit)/2):3020
%    addpoints(ellipse2, TransferOrbit(i, 1), TransferOrbit(i, 2), TransferOrbit(i, 3));
%    drawnow;
% end
% for i = 1:length(Pos)
%     addpoints(descent, Pos(i, 1), Pos(i, 2), Pos(i, 3));
%     drawnow;
% end


burnTime = Vehicle.mass/Vehicle.MassFlow * (1-exp(-abs(DeltaV)/Vehicle.Isp));    % s - duration of the deceleration burn
mass_afterBurn = Vehicle.mass - Vehicle.MassFlow * burnTime;                % kg - mass of Vehicle after burn

flight_path_angle = atan2(1 + e_t * sin(theta_entry), 1 + e_t * cos(theta_entry));
fprintf('flight path angle = %.4f°\n', rad2deg(flight_path_angle));
V_entry = V_entry * 1e3;   % we need the value in m/s
save('Data/orbit.mat', 'V_entry', 'flight_path_angle');



