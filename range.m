function rho = range(x, y, z, xs, ys, zs, theta)
%{
Function to calculate the range from a satellite to a ground station.
Ground station coordinates must be in ECEF [m] and satellite must be in ECI
[m].

Parameters
----------
x : double
    x-position of satellite in ECI [meters]
y : double
    y-position of satellite in ECI [meters]
z : double
    z-position of satellite in ECI [meters]
xs : double
    x-position of ground station in ECEF [meters]
ys : double
    y-position of ground station in ECEF [meters]
zs : double
    z-position of ground station in ECEF [meters]
theta : double
    greenwhich hour angle difference between ECI and ECEF frames

Returns
-------
rho : double
    range in meters
%}

rho = sqrt(x^2 + y^2 + z^2 + xs^2 + ys^2 + zs^2 - ...
    2*(x*xs + y*ys)*cos(theta) + ...
    2*(x*ys - y*xs)*sin(theta) - 2*z*zs);
end