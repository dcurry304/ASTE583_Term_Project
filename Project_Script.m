clc;clear;

% defining constants
constants = struct;
constants.Area = 3;                     % m^2
constants.Mass = 970;                   % kg
constants.rho0 = 3.614e-13;             % kg/m^3
constants.H = 88667.0;                  % m
constants.Re = 6378136.3;               % m
constants.r0 = 700000.0 + constants.Re; % m
constants.theta_dot = 7.2921158543e-5;  % rad/sec

%load observation data
obs = load('obs_data.mat');

%make the groundstation indexing easier
obs.groundstation(obs.groundstation==101,2) = 10;
obs.groundstation(obs.groundstation==337,2) = 13;
obs.groundstation(obs.groundstation==394,2) = 16;

% struct to save output data
out = struct; 

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
vals = [1e-6 1e-6 1e-6 ...
    1e-6 1e-6 1e-6 ...
    1e-20 1e-6 1e-6 ...
    1e10 1e10 1e10 ...
    1e-6 1e-6 1e-6 ...
    1e-6 1e-6 1e-6];
inv_P0_bar = diag(vals);

% weighing matrix
W = [1/(.01^2) 0; 0 1/(.001^2)];  % m and m/s  noise

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12); % seconds
theta = constants.theta_dot .* obs.time; % rads

iters = 4;

%-% BATCH PROCESSOR
j = 1;
while j < iters
    out.state_parameters(j,:) = x0';
    [t,sols] = ode45(@(t,Y) dynamics(Y, constants), ...
        obs.time, [x0; STM_IC(:)], options);
    x = sols(:,1);
    y = sols(:,2);
    z = sols(:,3);
    x_dot = sols(:,4);
    y_dot = sols(:,5);
    z_dot = sols(:,6);


    lambda = inv_P0_bar;
    N = inv_P0_bar * dx0_a_priori;

    for i = 1:length(obs.time)

        %Second column of grounstation is setup with correct index into x0
        idx = obs.groundstation(i,2);
        xs = x0(idx);
        ys = x0(idx+1);
        zs = x0(idx+2);

        H_tilde = H_tilde_xs1(theta(i),constants.theta_dot, ...
            x(i),x_dot(i),xs, ...
            y(i),y_dot(i),ys, ...
            z(i),z_dot(i),zs);
       
        % calculating range
        rho = range(x(i), y(i), z(i), ...
            xs, ys, zs, theta(i));
        
        % calculating range-rate
        rho_dot = range_rate(x(i), y(i), z(i), ...
            x_dot(i), y_dot(i), z_dot(i), ...
            xs, ys, zs, theta(i), constants.theta_dot);
        
        % calculating the observation residuals
        G = [rho; rho_dot];
        y_i = [obs.range(i); obs.range_rate(i)] - G;
        
        % getting the STM at timestep (i)
        % calculating the state-observation matrix and mapping
        % it to timestep (i) using the STM
        phi = reshape(sols(i, 19:end), 18, 18);
        H = H_tilde * phi;
        
        % updating normal equations
        lambda = lambda + (H.' * W * H);
        N = N + (H.' * W * y_i);
        
        % saving outputs for post-processing
        out.rho_residuals(i) = y_i(1);
        out.rho_dot_residuals(i) = y_i(2);
    end
    
    fprintf("rho rms = %f\n",rms(out.rho_residuals(:)))
    fprintf("rho dot rms = %f\n",rms(out.rho_dot_residuals(:)))

    % solving the normal equations
    R = chol(lambda);
    state_deviation  = R\(R'\N);

    % updating the initial state vector
    x0 = x0 + state_deviation;

    % shifting the a priori deviation vector by
    % the state deviation vector
    dx0_a_priori = dx0_a_priori - state_deviation;

    j = j+1;
end

