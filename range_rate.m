function rho_dot = range_rate(x, y, z, ...
    x_dot, y_dot, z_dot, ...
    xs, ys, zs, theta, theta_dot)

%{
Function to calcualte the range-rate from a satellite to a ground station.
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
x_dot : double
    x-velocity of satellite in ECI [meters/second]
y_dot : double
    y-velocity of satellite in ECI [meters/second]
z_dot : double
    z-velocity of satellite in ECI [meters/second]
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
rho_dot : double
    range-rate in meters/second
%}

rho = range(x, y, z, xs, ys, zs, theta);

rho_dot = (x*x_dot + y*y_dot + z*z_dot - ...
    (x_dot*xs + y_dot*ys)*cos(theta) + ...
    theta_dot*(x*xs + y*ys)*sin(theta) + ...
    (x_dot*ys - y_dot*xs)*sin(theta) + ...
    theta_dot*(x*ys - y*xs)*cos(theta) - ...
    z_dot*zs) / rho;
end