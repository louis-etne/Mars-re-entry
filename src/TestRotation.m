clc,clear, close all;
load('data/constants.mat', 'Mars')
%EARTH_SPHERE Generate an earth-sized sphere.
%   [X,Y,Z] = EARTH_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURFACE(X,Y,Z) produces a sphere equal to 
%   the radius of the earth in kilometers. The continents will be
%   displayed.
%
%   [X,Y,Z] = EARTH_SPHERE uses N = 50.
%
%   EARTH_SPHERE(N) and just EARTH_SPHERE graph the earth as a 
%   SURFACE and do not return anything.
%
%   EARTH_SPHERE(N,'mile') graphs the earth with miles as the unit rather
%   than kilometers. Other valid inputs are 'ft' 'm' 'nm' 'miles' and 'AU'
%   for feet, meters, nautical miles, miles, and astronomical units
%   respectively.
%
%   EARTH_SPHERE(AX,...) plots into AX instead of GCA.
% 
%  Examples: 
%    earth_sphere('nm') produces an earth-sized sphere in nautical miles
%
%    earth_sphere(10,'AU') produces 10 point mesh of the Earth in
%    astronomical units
%
%    h1 = gca;
%    earth_sphere(h1,'mile')
%    hold on
%    plot3(x,y,z)
%      produces the Earth in miles on axis h1 and plots a trajectory from
%      variables x, y, and z
%   Clay M. Thompson 4-24-1991, CBM 8-21-92.
%   Will Campbell, 3-30-2010
%   Copyright 1984-2010 The MathWorks, Inc. 
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
    
%     axlen = Mars.radius*1.2;
%     axis([-axlen, axlen, -axlen, axlen, -axlen, axlen])
    axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])

    xx = x; yy = y; zz = z;

    ax = gca;
    ax.Clipping = 'off';
    hold on
    % h = surf(X,Y,Z);
    % set(h,'EdgeColor','r')
    % set(h,'FaceCOlor', 'b');
    % alpha(h, 1)
%     t = 45;
%     st1 = sind( t*sind(t) );
%     ct1 = cosd( t*sind(t) );
    view(127.5,30)
    pause(1);
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
        if t ~= 0.1
            delete(h)
        end
    end




