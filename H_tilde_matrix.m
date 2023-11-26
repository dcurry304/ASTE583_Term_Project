function H_tilde = H_tilde_matrix(X,Xs,idx,theta,rho,rhodot,const)

%H_tilde = [d_rho/d_X; d_rhodot/d_X]
%H_tilde = 2x18 matrix

%reduce number of operations
cos_theta = cos(theta);
sin_theta = sin(theta);

H_tilde = zeros(2,const.sz);

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

%d_rhodot/d_x
H_tilde(2,1) = (X(4) + const.theta_dot*Xs(1)*sin_theta + const.theta_dot*Xs(2)*cos_theta)/rho ...
    - (rhodot*(X(1) - Xs(1)*cos_theta + Xs(2)*sin_theta))/rho^2;
 
%d_rhodot/d_y
H_tilde(2,2) = (X(5) + const.theta_dot*Xs(2)*sin_theta - const.theta_dot*Xs(1)*cos_theta)/rho ...
    - (rhodot*(X(2) - Xs(2)*cos_theta - Xs(1)*sin_theta))/rho^2;
 
%d_rhodot/d_z
H_tilde(2,3) = X(6)/rho - rhodot*(X(3)-Xs(3))/rho^2;
 
%d_rhodot/d_xdot
H_tilde(2,4) = (X(1) - Xs(1)*cos_theta + Xs(2)*sin_theta)/rho;
 
%d_rhodot/d_ydot
H_tilde(2,5) = (X(2) - Xs(2)*cos_theta - Xs(1)*sin_theta)/rho;
 
%d_rhodot/d_dot
H_tilde(2,6) = (X(3)-Xs(3))/rho;
 
%d_rhodot/d_xsi
H_tilde(2,idx) = (-X(4)*cos_theta + const.theta_dot*X(1)*sin_theta - X(5)*sin_theta - const.theta_dot*X(2)*cos_theta)/rho ...
    - rhodot*(Xs(1) - X(1)*cos_theta - X(2)*sin_theta)/rho^2;
 
%d_rhodot/d_ysi
H_tilde(2,idx+1) = (-X(5)*cos_theta + const.theta_dot*X(2)*sin_theta + X(4)*sin_theta + const.theta_dot*X(1)*cos_theta)/rho ...
    - rhodot*(Xs(2) - X(2)*cos_theta + X(1)*sin_theta)/rho^2;
 
%d_rhodot/d_zsi
H_tilde(2,idx+2) = -X(6)/rho - rhodot*(Xs(3)-X(3))/rho^2;
 
end