function [x, y, z] = SphereToCartesian(r, theta, phi)
    x = r*sin(phi)*cos(theta);
    y = r*sin(phi)*sin(theta);
    z = r*cos(phi);
end