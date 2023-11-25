function A = A_matrix(X, const)

r = norm(X(1:3));
v = X(4:6);
% velocity of the S/C wrt. the atmosphere
va = norm(v - cross(const.theta_dot*[0 0 1], X(1:3)).');

% atmospheric density calculation
rho = const.rho0 * exp(-(norm(r) - const.r0) / const.H);

drag_term = 1/2*const.CD*(const.Area/const.Mass)*rho;

A = zeros(18);

%d_xdot/d_xdot,same for y and z
A(1:3,4:6) = eye(3);

%d_xdotdot/d_x
A(4,1) = -const.mu/r^3*(1 - 3/2*const.J2*(const.Re/r)^2*(5*(X(3)/r)^2 -1)) ...
    + 3*const.mu*X(1)^2/r^5 *(1-5/2*const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-1)) ...
    + 1/2*const.CD*(const.Area/const.Mass)*rho*va*X(1)*(X(4)+const.theta_dot*X(2))/(r*const.H) ...
    -1/2*const.CD*(const.Area/const.Mass)*rho*(-const.theta_dot*X(5)+const.theta_dot^2*X(1))*(X(4)+const.theta_dot*X(2))/va;

%d_xdotdot/d_y
A(4,2) = 3*const.mu*X(1)*X(2)/r^5 *(1-5/2 *const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-1)) ...
    + drag_term*(va*X(2)*(X(4)+const.theta_dot*X(2))/(r*const.H) ...
    - (const.theta_dot*X(4)+const.theta_dot^2*X(2))*(X(4)+const.theta_dot*X(2))/va ...
    - va*const.theta_dot);

%d_xdotdot/d_z
A(4,3) = 3*const.mu*X(1)*X(3)/r^5 *(1-5/2 *const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-3)) ...
    + drag_term*va*X(3)*(X(4)+const.theta_dot*X(2))/(r*const.H);

%d)xdotdot/d_xdot
A(4,4) = -drag_term*((X(4)+const.theta_dot*X(2))^2/va ...
    -va);

%d_xdotdot/d_ydot
A(4,5) = -drag_term*(X(5)-const.theta_dot*X(1))*(X(4)+const.theta_dot*X(2))/va;

%d_xdotdot/d_zdot
A(4,6) = -drag_term*X(6)*(X(4)+const.theta_dot*X(2))/va;

%d_xdotdot/d_mu
A(4,7) = -X(1)/r^3*(1 - 3/2*const.J2*(const.Re/r)^2*(5*(X(3)/r)^2 -1));

%d_xdotxdot/d_J2
A(4,8) = 3/2 *const.mu*X(1)/r^3*((const.Re/r)^2*(5*(X(3)/r)^2 -1));

%d_xdotdot/d_CD
A(4,9) = -1/2*(const.Area/const.Mass)*rho*va*(X(4)+const.theta_dot*X(2));

%d_ydotdot/d_x
A(5,1) = 3*const.mu*X(1)*X(2)/r^5 *(1- 5/2*const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-1)) ...
    + drag_term*va*(X(5)-const.theta_dot*X(1))*X(1)/(r*const.H) ...
    -drag_term*(const.theta_dot^2*X(1)-const.theta_dot*X(5))*(X(5)-const.theta_dot*X(1))/va ...
    + drag_term*va*const.theta_dot;

%d_ydotdot/d_y
A(5,2) = -const.mu/r^3 *(1-3/2*const.J2*(const.Re/r)^2*(5*(X(3)/r)^2-1)) ...
    + 3*const.mu*X(2)^2/r^5 *(1-5/2*const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-1)) ...
    +drag_term*va*X(2)*(X(5)-const.theta_dot*X(1))/(r*const.H) ...
    -drag_term*(const.theta_dot*X(4)+const.theta_dot^2*X(2))*(X(5)-const.theta_dot*X(1))/va;

%d_ydotdot/d_z
A(5,3) = 3*const.mu*X(2)*X(3)/r^5 *(1- 5/2*const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-3)) ...
    + drag_term*va*X(3)*(X(5)-const.theta_dot*X(1))/(r*const.H);

%d_ydotdot/d_xdot
A(5,4) = -drag_term*(X(5)-const.theta_dot*X(1))*(X(4) + const.theta_dot*X(2))/va;

%d_ydotdot/d_ydot
A(5,5) = -drag_term*(X(5)-const.theta_dot*X(1))^2/va - drag_term*va;

%d_ydotdot/d_zdot
A(5,6) = -drag_term*X(6)*(X(4)-const.theta_dot*X(1))/va;

%d_ydotdot/d_mu
A(5,7) = -X(2)/r^3*(1-3/2*const.J2*(const.Re/r)^2*(5*(X(3)/r)^2-1));

%d_ydotdot/d_J2
A(5,8) = 3/2*const.mu*X(2)/r^3*(const.Re/r)^2*(5*(X(3)/r)^2-1);

%d_ydotdot/d_CD
A(5,9) = -1/2*(const.Area/const.Mass)*rho*va*(X(5)-const.theta_dot*X(1));

%d_zdotdot/d_x
A(6,1) = 3*const.mu*X(1)*X(3)/r^5 *(1 -5/2 *const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-3)) ...
    + drag_term*va*X(6)*X(1)/(r*const.H) - ...
    drag_term*X(6)*(const.theta_dot^2 *X(1) - const.theta_dot*X(5))/va;

%d_zdotdot/d_y
A(6,2) = 3*const.mu*X(2)*X(3)/r^5 *(1 -5/2 *const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-3)) ...
    + drag_term*va*X(6)*X(2)/(r*const.H) - ...
    drag_term*X(6)*(const.theta_dot^2 *X(1) + const.theta_dot^2*X(5))/va;

%d_zdotdot/d_z
A(6,3) = -const.mu/r^3*(1 - 3/2*const.J2*(const.Re/r)^2 *(5*(X(3)/r)^2 -3)) ...
    + 3*const.mu*X(3)^2/r^5 *(1 -5/2 *const.J2*(const.Re/r)^2*(7*(X(3)/r)^2-5)) ...
    + drag_term*va*X(3)*X(6)/(r*const.H);
    
%d_zdotdot/d_xdot
A(6,4) = -drag_term*X(6)*(X(4)+const.theta_dot*X(2))/va;

%d_zdotdot/d_ydot
A(6,5) = -drag_term*X(6)*(X(5)+const.theta_dot*X(1))/va;

%d_zdotdot/d_zdot
A(6,6) = -drag_term*X(6)^2/va - drag_term*va;

%d_zdotdot/d_mu
A(6,7) = -X(3)/r^3 *(1 - 3/2 *const.J2*(const.Re/r)^2*(5*(X(3)/r)^2 -3));

%d_zdotdot/d_J2
A(6,8) = 3/2*const.mu*X(3)/r^3 *(const.Re/r)^2*(5*(X(3)/r)^2 -3);

%d_zdotdot/d_CD
A(6,9) = -1/2*(const.Area/const.Mass)*rho*va*X(6);
end