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

%%%%%%%%%%%%%%%%%%%% initial conditions %%%%%%%%%%%%%%%%%%%%
mu = 3.986004415e14;                            % m^3/s^2
J2 = 1.082626925638815e-3;                      % N/A
CD = 2;                                         % N/A
ic.x0 = [r0; v0; mu; J2; CD; const.st1; const.st2; const.st3;reshape(eye(const.sz),const.sz*const.sz,1)];

% a-priori state deviation vector
ic.dx0_a_priori = zeros(const.sz,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a-priori Covariance
ic.inv_P0_bar = diag([1e-6*ones(6,1);1e-20;1e-6;1e-6;1e10*ones(3,1);1e-6*ones(6,1)]);

% weighing matrix
ic.W = [1/(0.01^2) 0; 0 1/(0.001^2)];  % m and m/s noise

fprintf("---------Range and range rate observations -----------\n");
%run case with both range and range rate observations
out = batch_processor(const,obs,ic,1);

% weighing matrix
ic.W = 1/(0.01^2);  % m noise

%run case with only range observations
fprintf("---------Range only observations -----------\n");
out2 = batch_processor(const,obs,ic,2);

% weighing matrix
ic.W = 1/(0.001^2);  % m/s noise

fprintf("---------Range-rate only observations -----------\n");
%run case with only range rate observations
out3 = batch_processor(const,obs,ic,3);