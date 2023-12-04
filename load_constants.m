function constants = load_constants()
% This function defines the constants assumed for this orbit determination
% problem.
% 
% Inputs
% ----------
% None
% 
% Outputs
% -------
% constants: structure of constants variables

constants = struct;
%state vector size (6 for s/c state + mu + J2 + CD + 3 per 3 stations)
constants.sz = 18;

%spacecraft cross sectional area
constants.Area = 3;                     % m^2
%spacecraft mass
constants.Mass = 970;                   % kg
%atmospheric density at intial time
constants.rho0 = 3.614e-13;             % kg/m^3
%scale height
constants.H = 88667.0;                  % m
%Earth's mean equatorial radius
constants.Re = 6378136.3;               % m
%spacecraft radius from center of earth at inital time
constants.r0 = 700000.0 + constants.Re; % m
%Earth's rotation rate about its axis
constants.theta_dot = 7.2921158543e-5;  % rad/sec

%initial station location (only station 1's location is exact, others are
%estimated by batch filter
constants.st1 = [-5127510.0; -3794160.0; 0.0];       % m ECEF
constants.st2 = [3860910.0; 3238490.0; 3898094.0];   % m ECEF
constants.st3 = [549505.0; -1380872.0; 6182197.0]; 

%options for the ode45 function
constants.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);% m ECEF

end