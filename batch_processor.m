function out = batch_processor(const,obs,ic,obs_opt)
% This function performs the batch filtering algorithm to determine the
% best estimate state at t0.
% 
% Inputs
% ----------
% const: structure of constants
% obs: structure of observation data
% ic: structure of initial conditions for this batch filter computation
% obs_opt: observation option to calculate:
%   1: use range and range-rate observations
%   2: use range only observations
%   3. use range-rate only observations
% 
% Outputs
% -------
% out: output structure with batch filter results
%  

%unload initial conditions structure to make it easy to access variables
x0 = ic.x0;
dx0_a_priori = ic.dx0_a_priori;
inv_P0_bar = ic.inv_P0_bar;
W = ic.W;

%-% BATCH PROCESSOR
j = 1;
iters = 5;
while j < iters
    %reset STM to Identity every pass
    x0(const.sz+1:end) = reshape(eye(const.sz),const.sz*const.sz,1);
    %solve for state vector using ode function
    [~,X] = ode45(@(t,Y)dynamics(t,Y, const), obs.time, x0, const.options);

    %use a priori estimate in normal equation
    M = inv_P0_bar;
    N = inv_P0_bar * dx0_a_priori;

    for i = 1:length(obs.time)
        %Second column of station is setup with correct index into x0
        idx = obs.station(i,2);
        %coordinates for station 1 is exact, others are estimated
        if idx == 10
            Xs = const.st1;
        else
            Xs = X(i,idx:idx+2);
        end
        
        % calculating range [m]
        rho = sqrt(norm(X(i,1:3))^2 + norm(Xs)^2 - ...
                2*(X(i,1)*Xs(1) + X(i,2)*Xs(2))*cos(obs.theta(i)) + ...
                2*(X(i,1)*Xs(2) - X(i,2)*Xs(1))*sin(obs.theta(i)) - 2*X(i,3)*Xs(3));
        
        % calculating range-rate [m/s]
        rho_dot = ( dot(X(i,1:3),X(i,4:6)) - ...
                  (X(i,4)*Xs(1) + X(i,5)*Xs(2))*cos(obs.theta(i)) + ...
                  const.theta_dot*(X(i,1)*Xs(1) + X(i,2)*Xs(2))*sin(obs.theta(i)) + ...
                  (X(i,4)*Xs(2) - X(i,5)*Xs(1))*sin(obs.theta(i)) + ...
                  const.theta_dot*(X(i,1)*Xs(2) - X(i,2)*Xs(1))*cos(obs.theta(i)) - ...
                  X(i,6)*Xs(3)) / rho;
        
        %Compute H_tilde matrix
        H_tilde = H_tilde_matrix(X(i,1:6),Xs,idx,obs.theta(i),rho,rho_dot,const,obs_opt);
        
        % calculating the observation residuals
        if obs_opt == 1
            G = [rho; rho_dot]; %[m,m/s]
            y_i = [obs.range(i); obs.range_rate(i)] - G; %[m,m/s]
        elseif obs_opt == 2
            G = rho;%[m]
            y_i = obs.range(i) - G;%[m]
        else
            G = rho_dot;%[m/s]
            y_i = obs.range_rate(i) - G;%[m/s]
        end
        
        % getting the STM at timestep (i)
        phi = reshape(X(i, const.sz+1:end), const.sz, const.sz);
        % calculating the state-observation matrix and mapping
        % it to timestep (i) using the STM
        H = H_tilde * phi;
        
        % updating normal equations
        M = M + (H.' * W * H);
        N = N + (H.' * W * y_i);
        
        % saving outputs for post-processing
        if obs_opt == 1 
            out.rho_residuals(j,i) = y_i(1); %[m]
            out.rho_dot_residuals(j,i) = y_i(2); %[m/s]
        elseif obs_opt == 2
            out.rho_residuals(j,i) = y_i(1);%[m]
        else
            out.rho_dot_residuals(j,i) = y_i(1);%[m/s]
        end
    end %time loop
    
    % solving the normal equations via Cholesky Decomposition
    [x_hat_0,P0] = Cholesky_Decomp(M,N);

    % updating the initial state vector for next iteration
    x0(1:const.sz) = x0(1:const.sz) + x_hat_0;

    % shifting the a priori deviation vector by
    % the state deviation vector for next iteration
    dx0_a_priori = dx0_a_priori - x_hat_0;
    
    %save outputs
    out.x_hat0(j,:) = x0(1:const.sz)';
    out.P0(j,:,:) = P0;
    out.traceP0(j) = trace(P0);
    out.sigmas(j,:) = diag(P0);
    
    %print results
    fprintf("Best estimate initial state:  r = [%.3e %.3e %.3e] m, v = [%.3e %.3e %.3e] m/s\n",out.x_hat0(j,1:6))
    fprintf("Best estimate initial sigmas: r = [%.3e %.3e %.3e] m, v = [%.3e %.3e %.3e] m/s\n",out.sigmas(j,1:6))
    if obs_opt == 1 || obs_opt == 2
        fprintf("rho rms = %f\n",rms(out.rho_residuals(j,:)))
    end
    if obs_opt == 1 || obs_opt == 3
        fprintf("rho dot rms = %f\n",rms(out.rho_dot_residuals(j,:)))
    end
    
    j = j+1;
end

%Propagate the best estimate initial state through to tf
x0(1:const.sz) = out.x_hat0(4,:)';
x0(const.sz+1:end) = reshape(eye(const.sz),const.sz*const.sz,1);
[~,X] = ode45(@(t,Y)dynamics(t,Y, const), obs.time, x0, const.options);
out.X = X(:,1:const.sz);
end