function Xdot = dynamics(~,X, const)
% This function is used to integrate both the state and the 
% state transition matrix of the two body problem.
% It integrates two equations:
% 
% 1) r_ddot = - mu / r^3 TODO: fix this 
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
r = X(1:3);
v = X(4:6);

mu = X(7);
J2 = X(8);
CD = X(9);

%get STM from state vec
phi = reshape(X(const.sz+1:end),const.sz,const.sz);

%compute A matrix
A = A_matrix(X,const);

% multiplying A matrix by the initial conditions
phi_dot = A * phi;

% velocity of the S/C wrt. the atmosphere
va = v - cross(const.theta_dot*[0 0 1], r).';

% atmospheric density calculation
rho = const.rho0 * exp(-(norm(r) - const.r0) / const.H);

% calculating drag acceleration
a_drag = (1/2) * CD * (const.Area/const.Mass) * rho * norm(va) .* va;

ax = -(mu / norm(r)^3)*r(1) ...
    - (((3*mu) / (2*norm(r)^5))*const.Re^2*J2*(1-(5*(r(3)/norm(r))^2))*r(1)) ...
    - a_drag(1);
ay = -(mu / norm(r)^3)*r(2) ...
    - (((3*mu) / (2*norm(r)^5))*const.Re^2*J2*(1-(5*(r(3)/norm(r))^2))*r(2)) ...
    - a_drag(2);
az = -(mu / norm(r)^3)*r(3) ...
    - (((3*mu) / (2*norm(r)^5))*const.Re^2*J2*(3-(5*(r(3)/norm(r))^2))*r(3)) ...
    - a_drag(3);

Xdot = [v; ax;ay;az; zeros(12,1); phi_dot(:)];
end
