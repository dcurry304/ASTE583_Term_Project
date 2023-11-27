function out = batch_processor(const,obs,ic,obs_opt)

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
        H_tilde = H_tilde_matrix(X(i,1:6),Xs,idx,obs.theta(i),rho,rho_dot,const,obs_opt);
        
        % calculating the observation residuals
        if obs_opt == 1
            G = [rho; rho_dot];
            y_i = [obs.range(i); obs.range_rate(i)] - G;
        elseif obs_opt == 2
            G = rho;
            y_i = obs.range(i) - G;
        else
            G = rho_dot;
            y_i = obs.range_rate(i) - G;
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
            out.rho_residuals(j,i) = y_i(1);
            out.rho_dot_residuals(j,i) = y_i(2);
        elseif obs_opt == 2
            out.rho_residuals(j,i) = y_i(1);
        else
            out.rho_dot_residuals(j,i) = y_i(1);
        end
    end
    % solving the normal equations via Cholesky Decomposition
    [x_hat_0,P0] = Cholesky_Decomp(M,N);

    % updating the initial state vector
    x0(1:const.sz) = x0(1:const.sz) + x_hat_0;

    % shifting the a priori deviation vector by
    % the state deviation vector
    dx0_a_priori = dx0_a_priori - x_hat_0;
    
    %save outputs
    out.x_hat0(j,:) = x0(1:const.sz)';
    out.P0(j,:,:) = P0;
    out.sigmas(j,:) = diag(P0);
    
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
end