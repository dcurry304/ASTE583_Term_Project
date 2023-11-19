function constants = load_constants()

% defining constants
constants = struct;
%state vector size (6 for s/c state + mu + J2 + CD + 3 per 3 stations)
constants.sz = 18;

constants.Area = 3;                     % m^2
constants.Mass = 970;                   % kg
constants.rho0 = 3.614e-13;             % kg/m^3
constants.H = 88667.0;                  % m
constants.Re = 6378136.3;               % m
constants.r0 = 700000.0 + constants.Re; % m
constants.theta_dot = 7.2921158543e-5;  % rad/sec

constants.mu = 3.986004415e14;                            % m^3/s^2
constants.J2 = 1.082626925638815e-3;                      % N/A
constants.CD = 2;                                         % N/A
constants.st1 = [-5127510.0; -3794160.0; 0.0];       % m ECEF
constants.st2 = [3860910.0; 3238490.0; 3898094.0];   % m ECEF
constants.st3 = [549505.0; -1380872.0; 6182197.0]; 

constants.options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);% m ECEF

end