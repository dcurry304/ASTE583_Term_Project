function phi = dynamics(X, const)
% This function is used to integrate both the state and the 
% state transition matrix of the two body problem.
% It integrates two equations:
% 
% 1) r_ddot = - mu / r^3
% 2) phi_dot = dx_dot/ dx * phi
% 
% Parameters
% ----------
% mu
% X
% 
% Returns
% -------
% phi
% 

% getting the state vector -> [x, y, z, vx, vy, vz] 
r = X(1:3);
v = X(4:6);

x = r(1);
y = r(2);
z = r(3);

x_dot = v(1);
y_dot = v(2);
z_dot = v(3);

mu = X(7);
J2 = X(8);
CD = X(9);

Area = const.Area;
H = const.H;
Re = const.Re;
r0 = const.r0;
m = const.Mass;
theta_dot = const.theta_dot;
rho0 = const.rho0;

% getting the stm (phi)
stm_ic = X(19:end);

A = A_Matrix(Area,CD,H,J2,Re,m,mu,r0,rho0,theta_dot,x,x_dot,y,y_dot,z,z_dot);

% multiplying A matrix by the initial conditions
stm_dt = A * reshape(stm_ic,18,18);

% % % % % % % % % % % % % % Drag Acceleration % % % % % % % % % % % % % %
% velocity of the S/C wrt. the atmosphere
va = v - cross(theta_dot*[0 0 1], r).';

% atmospheric density calculation
rho = rho0 * exp(-(norm(r) - r0) / H);

% calculating drag acceleration
a_drag = (1/2) * CD * (Area/m) * rho * norm(va) .* va;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

ax = -(mu / norm(r)^3)*x ...
    - (((3*mu) / (2*norm(r)^5))*Re^2*J2*(1-(5*(z/norm(r))^2))*x) ...
    - a_drag(1);
ay = -(mu / norm(r)^3)*y ...
    - (((3*mu) / (2*norm(r)^5))*Re^2*J2*(1-(5*(z/norm(r))^2))*y) ...
    - a_drag(2);
az = -(mu / norm(r)^3)*z ...
    - (((3*mu) / (2*norm(r)^5))*Re^2*J2*(3-(5*(z/norm(r))^2))*z) ...
    - a_drag(3);

phi = [v; ax; ay; az; 0; 0; 0; [0; 0; 0]; [0; 0; 0]; [0; 0; 0]; stm_dt(:)];
end
