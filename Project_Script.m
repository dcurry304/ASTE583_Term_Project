%ASTE 583 Final Project: Batch processor for the 2 body problem with drag
%and J2 effects. This driver script 
%   1. loads the inputs (constants structure parameters, batch filter parameters, and observation data)
%   2. calls the batch processor 3 times (range + range-rate, range only, range-rate only
%observations)
%   3. prints out the results of each batch filter call
%   4. calls the Kalman filter function once and plots its results
%
close all;clear;

%load in constants 
const = load_constants();

%initial spacecraft position and velocity in the inertial coordinate system
r0 = [757700.0; 5222607.0; 4851500.0]; % [m]
v0 = [2213.21; 4678.34; -5371.30];% [m/s]

%load observation data
load('obs_data.mat');
%the angle between the inertial frame and body-fixed frame.
obs.theta = const.theta_dot .* obs.time; % rads

%index of the station i in the state vector
obs.station(obs.station==101,2) = 10;
obs.station(obs.station==337,2) = 13;
obs.station(obs.station==394,2) = 16;

%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%
%earth gravitational parameter, to be estimated by filter
mu = 3.986004415e14;                            % m^3/s^2
%earth J2 parameter, to be estimated by filter
J2 = 1.082626925638815e-3;                      % N/A
%spacecraft coefficient of drag, to be estimated by filter
CD = 2;                                         % N/A
%initial state vector
ic.x0 = [r0; v0; mu; J2; CD; const.st1; const.st2; const.st3;reshape(eye(const.sz),const.sz*const.sz,1)];

% a-priori state deviation vector
ic.dx0_a_priori = zeros(const.sz,1);
% a-priori Covariance
ic.inv_P0_bar = diag([1e-6*ones(6,1);1e-20;1e-6;1e-6;1e10*ones(3,1);1e-6*ones(6,1)]);
ic.P0_bar = diag([1e6*ones(6,1);1e20;1e6;1e6;1e-10*ones(3,1);1e6*ones(6,1)]);

% weighing matrix
ic.W = [1/(0.01^2) 0; 0 1/(0.001^2)];  % m and m/s noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf("---------Range and range rate observations -----------\n");
%run case with both range and range rate observations
out = batch_processor(const,obs,ic,1);

%run case with only range observations
fprintf("---------Range only observations -----------\n");
% weighing matrix
ic.W = 1/(0.01^2);  % m noise

out2 = batch_processor(const,obs,ic,2);

fprintf("---------Range-rate only observations -----------\n");
% weighing matrix
ic.W = 1/(0.001^2);  % m/s noise

%run case with only range rate observations
out3 = batch_processor(const,obs,ic,3);

%Plot results
for i = 1:length(out.rho_dot_residuals(:,1))
    figure; 
    sgtitle(sprintf('Residuals for iteration %d',i));
    
    subplot(2,1,1); hold on;
    xlabel('Time (s)'); ylabel('Range Residuals (m)');
    plot(obs.time,out.rho_residuals(i,:))
    plot(obs.time,out2.rho_residuals(i,:))
    legend('Range + Range-Rate obs','Range only obs');
    
    subplot(2,1,2); hold on;
    xlabel('Time (s)'); ylabel('Range-rate Residuals (m)');
    plot(obs.time,out.rho_dot_residuals(i,:))
    plot(obs.time,out3.rho_dot_residuals(i,:))
    legend('Range + Range-Rate obs','Range-Rate only obs');
end

% weighing matrix
ic.W = [1/(0.01^2) 0; 0 1/(0.001^2)];  % m and m/s noise
%extra orbit determination task
fprintf("Kalman filter!\n");
outK = Kalman_filter(const,obs,ic,1);

figure; 
sgtitle('Kalman Filter Position deviation');
subplot(3,1,1); hold on;
xlabel('Time (s)'); ylabel('\deltaX (m)');
subplot(3,1,2); hold on;
xlabel('Time (s)'); ylabel('\deltaY (m)');
subplot(3,1,3); hold on;
xlabel('Time (s)'); ylabel('\deltaZ (m)');
for i=1:3
    subplot(3,1,i); hold on;
    plot(obs.time,outK.x_hat(:,i)-out.X(:,i))
    plot(obs.time,outK.dx_hat(:,i))
    legend('Batch', 'Reference');
end

figure; 
sgtitle('Kalman Filter Velocity deviation');
subplot(3,1,1); hold on;
xlabel('Time (s)'); ylabel('\deltaXdot (m)');
subplot(3,1,2); hold on;
xlabel('Time (s)'); ylabel('\deltaYdot (m)');
subplot(3,1,3); hold on;
xlabel('Time (s)'); ylabel('\deltaZdot (m)');
for i=4:6
    subplot(3,1,i-3); hold on;
    plot(obs.time,outK.x_hat(:,i)-out.X(:,i))
    plot(obs.time,outK.dx_hat(:,i))
    legend('Batch', 'Reference');
end
