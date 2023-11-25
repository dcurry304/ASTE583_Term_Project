function [x_hat, P] = Cholesky_Decomp(M,N)
%This function compute the best estimate state deviation, x_hat, and
%associated covariance matrix, P, using the Cholesky Decomposition method.
%The normal equation for the least squares estimation problem with a priori
%info is given below:
%   (H^T*W*H + inv_P_bar)*x_hat = H^T*W*y + inv_P_bar*x_bar
%which can be expressed as:
%   M*xhat = N
%
%Steps:
%   1. Solve for R in the equation R^T*R = M
%   2. Compute z given R^T*z = N
%   3. Compute x_hat given R*x^hat = z
%   4. Compute S = R^-1
%   5. Compute P = S*S^T

%Step 1: solve for R
R = zeros(size(M));
for i = 1:length(M)
    summ = 0;
    for k = 1:i-1
        summ = summ + R(k,i)^2;
    end
    R(i,i) = sqrt(M(i,i) - summ); 
    
    for j = i+1:length(M) 
        summ = 0;
        for k = 1:i-1
            summ = summ + R(k,i)*R(k,j);
        end
        R(i,j) = (M(i,j) - summ)/R(i,i); 
    end
end

%Step 2: solve for z
z = zeros(length(M),1);
for i = 1:length(M)
    summ = 0;
    for j = 1:i-1
        summ = summ + R(j,i)*z(j);
    end
   z(i,1) = (N(i)-summ)/R(i,i);
end

%Step 3: Compute x_hat
x_hat = zeros(length(M),1);
for i =length(M):-1:1
    summ = 0;
    for j = i+1:length(M)
        summ = summ + R(i,j)*x_hat(j);
    end
    x_hat(i) = (z(i)- summ)/R(i,i);
end

%Step 4: Compute S = R^-1
S = zeros(size(M));
for i = 1:length(M)
    S(i,i) = 1/R(i,i);
end
for i = 1:length(M) 
    for j = i+1:length(M)
        summ = 0;
        for k = i:j-1
            summ = summ + R(k,j)*S(i,k);
        end
        S(i,j) = -S(j,j)*summ;
    end
end

%Step 5: Compute P = S*S^T
P = S*S';

end