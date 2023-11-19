function Xdot = dynamics(X, const)
% This function is used to integrate both the state and the 
% state transition matrix of the two body problem.
% It integrates two equations:
% 
% 1) r_ddot = - const.mu / r^3
% 2) phi_dot = dx_dot/ dx * phi
% 
% Parameters
% ----------
% const.mu
% X
% 
% Returns
% -------
% phi
%  
r = X(1:3);
v = X(4:6);

%get STM from state vec
phi = reshape(X(const.sz+1:end),const.sz,const.sz);

A = A_Matrix(const.Area,const.CD,const.H,const.J2,const.Re,const.Mass,const.mu,const.r0,const.rho0,const.theta_dot,r(1),v(1),r(2),v(2),r(3),v(3));

% multiplying A matrix by the initial conditions
phi_dot = A * phi;

% % % % % % % % % % % % % % Drag Acceleration % % % % % % % % % % % % % %
% velocity of the S/C wrt. the atmosphere
va = v - cross(const.theta_dot*[0 0 1], r).';

% atmospheric density calculation
rho = const.rho0 * exp(-(norm(r) - const.r0) / const.H);

% calculating drag acceleration
a_drag = (1/2) * const.CD * (const.Area/const.Mass) * rho * norm(va) .* va;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

ax = -(const.mu / norm(r)^3)*r(1) ...
    - (((3*const.mu) / (2*norm(r)^5))*const.Re^2*const.J2*(1-(5*(r(3)/norm(r))^2))*r(1)) ...
    - a_drag(1);
ay = -(const.mu / norm(r)^3)*r(2) ...
    - (((3*const.mu) / (2*norm(r)^5))*const.Re^2*const.J2*(1-(5*(r(3)/norm(r))^2))*r(2)) ...
    - a_drag(2);
az = -(const.mu / norm(r)^3)*r(3) ...
    - (((3*const.mu) / (2*norm(r)^5))*const.Re^2*const.J2*(3-(5*(r(3)/norm(r))^2))*r(3)) ...
    - a_drag(3);

Xdot = [v; ax; ay; az; zeros(12,1); phi_dot(:)];
end
