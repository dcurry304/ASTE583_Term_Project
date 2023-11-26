clc;clear;

%load in constants 
const = load_constants();

r0 = [757700.0; 5222607.0; 4851500.0];           % m ECI
v0 = [2213.21; 4678.34; -5371.30];               % m/s ECI

%load observation data
load('obs_data.mat');
obs.theta = const.theta_dot .* obs.time; % rads

%make the station indexing easier
obs.station(obs.station==101,2) = 10;
obs.station(obs.station==337,2) = 13;
obs.station(obs.station==394,2) = 16;

% struct to save output data
out = struct; 

%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%
% state vector
STM_IC = eye(const.sz);
STM_IC = reshape(STM_IC,const.sz*const.sz,1);
x0 = [r0; v0; const.mu; const.J2; const.CD; const.st1; const.st2; const.st3;STM_IC];

% a-priori state deviation vector
dx0_a_priori = zeros(const.sz,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a-priori Covariance
inv_P0_bar = diag([1e-6*ones(6,1);1e-20;1e-6;1e-6;1e10*ones(3,1);1e-6*ones(6,1)]);

% weighing matrix
W = [1/(0.01^2) 0; 0 1/(0.001^2)];  % m and m/s noise

%-% BATCH PROCESSOR
j = 1;
iters = 5;
while j < iters
    out.state(j,:) = x0(1:const.sz)';
    %reset STM to Identity every pass
    x0(const.sz+1:end) = STM_IC;
    %solve for state vector using ode function
    [t,X] = ode45(@(t,Y) dynamics(Y, const), obs.time, x0, const.options);

    M = inv_P0_bar;
    N = inv_P0_bar * dx0_a_priori;

    for i = 1:length(obs.time)
        %Second column of station is setup with correct index into x0
        idx = obs.station(i,2);
        Xs = x0(idx:idx+2);
        
        % calculating range
        rho = sqrt(X(i,1)^2 + X(i,2)^2 + X(i,3)^2 + Xs(1)^2 + Xs(2)^2 + Xs(3)^2 - ...
                2*(X(i,1)*Xs(1) + X(i,2)*Xs(2))*cos(obs.theta(i)) + ...
                2*(X(i,1)*Xs(2) - X(i,2)*Xs(1))*sin(obs.theta(i)) - 2*X(i,3)*Xs(3));
        
        % calculating range-rate
        rho_dot = (X(i,1)*X(i,4) + X(i,2)*X(i,5) + X(i,3)*X(i,6) - ...
                  (X(i,4)*Xs(1) + X(i,5)*Xs(2))*cos(obs.theta(i)) + ...
                  const.theta_dot*(X(i,1)*Xs(1) + X(i,2)*Xs(2))*sin(obs.theta(i)) + ...
                  (X(i,4)*Xs(2) - X(i,5)*Xs(1))*sin(obs.theta(i)) + ...
                  const.theta_dot*(X(i,1)*Xs(2) - X(i,2)*Xs(1))*cos(obs.theta(i)) - ...
                  X(i,6)*Xs(3)) / rho;
        
        %Compute H_tilde matrix
        H_tilde = H_tilde_matrix(X(i,1:6),Xs,idx,obs.theta(i),rho,rho_dot,const);
        
        % calculating the observation residuals
        G = [rho; rho_dot];
        y_i = [obs.range(i); obs.range_rate(i)] - G;
        
        % getting the STM at timestep (i)
        phi = reshape(X(i, const.sz+1:end), const.sz, const.sz);
        % calculating the state-observation matrix and mapping
        % it to timestep (i) using the STM
        H = H_tilde * phi;
        
        % updating normal equations
        M = M + (H.' * W * H);
        N = N + (H.' * W * y_i);
        
        % saving outputs for post-processing
        out.rho_residuals(i) = y_i(1);
        out.rho_dot_residuals(i) = y_i(2);
    end
    
    fprintf("rho rms = %f\n",rms(out.rho_residuals(:)))
    fprintf("rho dot rms = %f\n",rms(out.rho_dot_residuals(:)))

    %R = chol(M);
    %state_deviation  = R\(R'\N);
    % solving the normal equations via Cholesky Decomposition
    [x_hat_0,P] = Cholesky_Decomp(M,N);

    % updating the initial state vector
    x0(1:const.sz) = x0(1:const.sz) + x_hat_0;

    % shifting the a priori deviation vector by
    % the state deviation vector
    dx0_a_priori = dx0_a_priori - x_hat_0;

    j = j+1;
end

