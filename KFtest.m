clc;clear;

%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%
% state vector
STM_IC = eye(18);
x0 = [3;2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a-priori Covariance
P_bar = [1 0; 0 1];

% weighing matrix
R = [2 0; 0 3/4];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
time_span = obs.time;                           % seconds
theta = constants.theta_dot .* time_span;       % rads

%% KALMAN FILTER
P = P_bar;
i = 1;
iters = 2;
while i < iters
    
    
    fprintf("%i\n", i)
    
    phi = [1 1; 0 1];

    H_tilde = [0 1; 0.5 0.5];
    

    % time update
    x_bar = phi*x0;
    P_bar = phi*P*phi.';
    
    y_i = [6; 4];


    K = P_bar*H_tilde.' * inv(H_tilde * P_bar * H_tilde.' + R);

    % measurement update
    x0 = x_bar + K*(y_i-(H_tilde*x_bar));
    P = (eye(2) - K*H_tilde)*P_bar;

    i = i + 1;
end