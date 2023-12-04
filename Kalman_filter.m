function out = Kalman_filter(const,obs,ic,obs_opt)

%unload initial conditions structure to make it easy to access variables
dx_hat = ic.dx0_a_priori;
out.x_hat(1,:) =ic.x0(1:const.sz);
R = inv(ic.W);
P = ic.P0_bar;

%-% Sequential Processor
%solve for state vector using ode function
[~,X] = ode45(@(t,Y)dynamics(t,Y, const), obs.time, ic.x0, const.options);
for i = 2:length(obs.time)

    %Second column of station is setup with correct index into x0
    idx = obs.station(i,2);
    %coordinates for station 1 is exact, others are estimated
    if idx == 10
        Xs = const.st1;
    else
        Xs = X(i,idx:idx+2);
    end

    % calculating range
    rho = sqrt(norm(X(i,1:3))^2 + norm(Xs)^2 - ...
            2*(X(i,1)*Xs(1) + X(i,2)*Xs(2))*cos(obs.theta(i)) + ...
            2*(X(i,1)*Xs(2) - X(i,2)*Xs(1))*sin(obs.theta(i)) - 2*X(i,3)*Xs(3));

    % calculating range-rate
    rho_dot = ( dot(X(i,1:3),X(i,4:6)) - ...
              (X(i,4)*Xs(1) + X(i,5)*Xs(2))*cos(obs.theta(i)) + ...
              const.theta_dot*(X(i,1)*Xs(1) + X(i,2)*Xs(2))*sin(obs.theta(i)) + ...
              (X(i,4)*Xs(2) - X(i,5)*Xs(1))*sin(obs.theta(i)) + ...
              const.theta_dot*(X(i,1)*Xs(2) - X(i,2)*Xs(1))*cos(obs.theta(i)) - ...
              X(i,6)*Xs(3)) / rho;
          
    % getting the STM at timestep (i)
    phi = reshape(X(i, const.sz+1:end), const.sz, const.sz);
    
    % time update
    x_bar = phi*dx_hat(1:const.sz);
    P_bar = phi*P*phi.';
    
    % calculating the observation residuals
    if obs_opt == 1
        G = [rho; rho_dot];
        y_i = [obs.range(i); obs.range_rate(i)] - G;
        out.rho_residuals(i) = y_i(1);
        out.rho_dot_residuals(i) = y_i(2);
    elseif obs_opt == 2
        G = rho;
        y_i = obs.range(i) - G;
        out.rho_residuals(i) = y_i(1);
    else
        G = rho_dot;
        y_i = obs.range_rate(i) - G;
        out.rho_dot_residuals(i) = y_i(1);
    end
    
    %Compute H_tilde matrix
    H_tilde = H_tilde_matrix(X(i,1:6),Xs,idx,obs.theta(i),rho,rho_dot,const,obs_opt);

    %Compute Kalman Gain
    K = P_bar*H_tilde.' * inv(H_tilde * P_bar * H_tilde.' + R);
    
    % measurement update
    dx_hat = x_bar + K*(y_i-(H_tilde*x_bar));
    P = (eye(18) - K*H_tilde)*P_bar;
    
    %save outputs
    out.dx_hat(i,:) = dx_hat(1:const.sz)';
    out.x_hat(i,:) = X(i,1:const.sz)' + dx_hat(1:const.sz);
    
    out.P(i,:,:) = P;
    out.traceP(i) = trace(P);
    out.sigmas(i,:) = diag(P);
    
end
out.X = X(:,1:const.sz);

% fprintf("Best estimate initial state:  r = [%.3e %.3e %.3e] m, v = [%.3e %.3e %.3e] m/s\n",out.x_hat0(j,1:6))
% fprintf("Best estimate initial sigmas: r = [%.3e %.3e %.3e] m, v = [%.3e %.3e %.3e] m/s\n",out.sigmas(j,1:6))
if obs_opt == 1 || obs_opt == 2
    fprintf("rho rms = %f\n",rms(out.rho_residuals))
end
if obs_opt == 1 || obs_opt == 3
    fprintf("rho dot rms = %f\n",rms(out.rho_dot_residuals))
end
