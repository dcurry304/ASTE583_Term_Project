function H_tilde = H_tilde_matrix(X,Xs,idx,theta,rho,rhodot,const,obs_opt)  
% This function calculates the H_tilde matrix for the 2 body problem with
% drag and J2 effects. 
% 
% Inputs
% ----------
% X: 1x6 state vector, here only r and v are needed from state
%   [x y z xdot ydot zdot]^T
% Xs: station location for this timestep's observations
% idx: index into the state vector for this station. Station 1, idx=10,
%   Station 2, index = 13, Station 3, index = 16
% rho: calculated range from the station to the spacecraft (not observed
%   range)
% rhodot: calculated range-rate from the station to the spacecraft (not
%   observed range-rate)
% const: structure of constants
% obs_opt: observation option to calculate:
%   1: use range and range-rate observations
%   2: use range only observations
%   3. use range-rate only observations
% 
% Outputs
% -------
% H_tilde: 2x18 observation-state matrix, filled with zeros except for
% specified members defined below
%   H_tilde = [d_rho/d_X; d_rhodot/d_X]
%  

%determine which row to place rhodot partials, since obs_opt ==3 does not
%use rho observations
rhodot_ind = 2;
if obs_opt == 3
    H_tilde = zeros(1,const.sz);
    rhodot_ind = 1;
elseif obs_opt == 2
    H_tilde = zeros(1,const.sz); 
else 
    H_tilde = zeros(2,const.sz);
end

%reduce number of operations
cos_theta = cos(theta);
sin_theta = sin(theta);

if obs_opt == 1 || obs_opt == 2
    %d_rho/d_x
    H_tilde(1,1) = (X(1) - Xs(1)*cos_theta + Xs(2)*sin_theta)/rho;

    %d_rho/d_y
    H_tilde(1,2) = (X(2) - Xs(2)*cos_theta - Xs(1)*sin_theta)/rho;

    %d_rho/d_z
    H_tilde(1,3) = (X(3)-Xs(3))/rho;

    %d_rho/d_xsi
    H_tilde(1,idx) = (Xs(1)- X(1)*cos_theta - X(2)*sin_theta)/rho;

    %d_rho/d_ysi
    H_tilde(1,idx+1) = (Xs(2)-X(2)*cos_theta + X(1)*sin_theta)/rho;

    %d_rho/d_zsi
    H_tilde(1,idx+2) = (Xs(3)-X(3))/rho;
end
if obs_opt == 1 || obs_opt == 3
    %d_rhodot/d_x
    H_tilde(rhodot_ind,1) = (X(4) + const.theta_dot*Xs(1)*sin_theta + ...
        const.theta_dot*Xs(2)*cos_theta)/rho ...
        - (rhodot*(X(1) - Xs(1)*cos_theta + Xs(2)*sin_theta))/rho^2;

    %d_rhodot/d_y
    H_tilde(rhodot_ind,2) = (X(5) + const.theta_dot*Xs(2)*sin_theta - ...
        const.theta_dot*Xs(1)*cos_theta)/rho ...
        - (rhodot*(X(2) - Xs(2)*cos_theta - Xs(1)*sin_theta))/rho^2;

    %d_rhodot/d_z
    H_tilde(rhodot_ind,3) = X(6)/rho - rhodot*(X(3)-Xs(3))/rho^2;

    %d_rhodot/d_xdot
    H_tilde(rhodot_ind,4) = (X(1) - Xs(1)*cos_theta + Xs(2)*sin_theta)/rho;

    %d_rhodot/d_ydot
    H_tilde(rhodot_ind,5) = (X(2) - Xs(2)*cos_theta - Xs(1)*sin_theta)/rho;

    %d_rhodot/d_dot
    H_tilde(rhodot_ind,6) = (X(3)-Xs(3))/rho;

    %d_rhodot/d_xsi
    H_tilde(rhodot_ind,idx) = (-X(4)*cos_theta + const.theta_dot*X(1)*sin_theta ...
        - X(5)*sin_theta - const.theta_dot*X(2)*cos_theta)/rho ...
        - rhodot*(Xs(1) - X(1)*cos_theta - X(2)*sin_theta)/rho^2;

    %d_rhodot/d_ysi
    H_tilde(rhodot_ind,idx+1) = (-X(5)*cos_theta + const.theta_dot*X(2)*sin_theta ...
        + X(4)*sin_theta + const.theta_dot*X(1)*cos_theta)/rho ...
        - rhodot*(Xs(2) - X(2)*cos_theta + X(1)*sin_theta)/rho^2;

    %d_rhodot/d_zsi
    H_tilde(rhodot_ind,idx+2) = -X(6)/rho - rhodot*(Xs(3)-X(3))/rho^2;
end
 
end