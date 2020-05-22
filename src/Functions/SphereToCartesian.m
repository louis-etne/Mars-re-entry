function [Pos] = SphereToCartesian(r, theta, phi)
    x = r.*sind(phi).*cosd(theta);
    y = r.*sind(phi).*sind(theta);
    z = r.*cosd(phi);
    Pos(:, 1) = x;
    Pos(:, 2) = y;
    Pos(:, 3) = z;
end