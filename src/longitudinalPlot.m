longitude_init = 0;
latitude_init = 0;
downrange = linspace(0, 0.05, length(h.data(:,1)));
downrange = rad2deg(downrange);
x = zeros(1, length(downrange));

plot3(downrange, x, h.data(:, 1))
title('Position')
xlabel('Longitude (°)');
ylabel('Latitude (°)');
zlabel('Altitude (m)')