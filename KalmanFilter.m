

% defining constants
constants = struct;
constants.Area = 3;                             % m^2
constants.Mass = 970;                           % kg
constants.rho0 = 3.614e-13;                     % kg/m^3
constants.H = 88667.0;                          % m
constants.Re = 6378136.3;                       % m
constants.r0 = 700000.0 + constants.Re;         % m
constants.theta_dot = 7.2921158543e-5;          % rad/sec

obs = load("data.mat");

mu = 3.986004415e14;                            % m^3/s^2
J2 = 1.082626925638815e-3;                      % N/A
CD = 2;                                         % N/A
station1 = [-5127510.0; -3794160.0; 0.0];       % m ECEF
station2 = [3860910.0; 3238490.0; 3898094.0];   % m ECEF
station3 = [549505.0; -1380872.0; 6182197.0];   % m ECEF
r = [757700.0; 5222607.0; 4851500.0];           % m ECI
v = [2213.21; 4678.34; -5371.30];               % m/s ECI

%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%
% state vector
STM_IC = eye(18);
x0 = [r; v; mu; J2; CD; station1; station2; station3];
% a-priori state deviation vector
dx0_a_priori = zeros(18,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a-priori Covariance
vals = [1e6 1e6 1e6 ...
    1e6 1e6 1e6 ...
    1e20 1e6 1e6 ...
    1e-10 1e-10 1e-10 ...
    1e6 1e6 1e6 ...
    1e6 1e6 1e6];
P_bar = diag(vals);

% weighing matrix
W = [1/(.01^2) 0; 0 1/(.001^2)];                % m and m/s  noise
R = inv(W);

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
time_span = obs.time;                           % seconds
theta = constants.theta_dot .* time_span;       % rads

%% KALMAN FILTER
P = P_bar;
i = 1;
iters = length(obs.time)-1;
while i < iters
    
    fprintf("%i\n", i)
    [t,sols] = ode45(@(t,Y) dynamics(Y, constants), ...
        [time_span(i) (time_span(i)+time_span(i+1))/2 time_span(i+1)], [x0; STM_IC(:)], options);
    
    x = sols(1,1);
    y = sols(1,2);
    z = sols(1,3);
    x_dot = sols(1,4);
    y_dot = sols(1,5);
    z_dot = sols(1,6);
    % mapping the STM to the correct time
    % from ode we get three timesteps t0 t_mid ti+1
    % this multiplies STM of the last two timesteps to map it to the first
    phi = reshape(sols(end-1, 19:end), 18, 18)*reshape(sols(end, 19:end), 18, 18);

    ground_station =  obs.groundstation(i);

    if ground_station == 101
        xs = x0(10);
        ys = x0(11);
        zs = x0(12);

        H_tilde = H_tilde_xs1(theta(i),constants.theta_dot, ...
            x,x_dot,xs, ...
            y,y_dot,ys, ...
            z,z_dot,zs);
    end

    if ground_station == 337
        xs = x0(13);
        ys = x0(14);
        zs = x0(15);

        H_tilde = H_tilde_xs2(theta(i),constants.theta_dot, ...
            x,x_dot,xs, ...
            y,y_dot,ys, ...
            z,z_dot,zs);
    end

    if ground_station == 394
        xs = x0(16);
        ys = x0(17);
        zs = x0(18);

        H_tilde = H_tilde_xs3(theta(i),constants.theta_dot, ...
            x,x_dot,xs, ...
            y,y_dot,ys, ...
            z,z_dot,zs);
    end
    

    % time update
    x_bar = phi*x0;
    P_bar = phi*P*phi.';
    
    % calculating range
    rho = range(x, y, z, ...
        xs, ys, zs, theta(i));
    
    % calculating range-rate
    rho_dot = range_rate(x, y, z, ...
        x_dot, y_dot, z_dot, ...
        xs, ys, zs, theta(i), constants.theta_dot);
    
    % calculating the observation residuals
    G = [rho; rho_dot];
    y_i = [obs.range(i); obs.range_rate(i)] - G;


    K = P_bar*H_tilde.' * inv(H_tilde * P_bar * H_tilde.' + R);

    % measurement update
    x0 = x_bar + K*(y_i-(H_tilde*x_bar));
    P = (eye(18) - K*H_tilde)*P_bar;

    i = i + 1;
end

function Y = inv_chol(L)
% Matrix Inversion using Cholesky Decomposition

N = size(L, 1) ;
Y = zeros(N, N) ;
% Work with the upper triangular matrix
R = L' ;
% Construct the auxillary diagonal matrix S = 1/rii
S = inv(diag(diag(R))) ;

for j=N:-1:1
    for i=j:-1:1
        Y(i,j) = S(i,j) - R(i,i+1:end)*Y(i+1:end,j) ;
        Y(i,j) = Y(i,j)/R(i,i) ;
        % Write out the symmetric element
        Y(j,i) = conj(Y(i,j)) ;
    end
end
end

%% HELPER FUNCTIONS
function rho = range(x, y, z, xs, ys, zs, theta)
%{
Function to calcualte the range from a satellite to a ground station.
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

A = A_Matrix(Area,CD,H,J2,Re,m,mu,r0,rho0,theta_dot, ...
    x,x_dot,y,y_dot,z,z_dot);

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
