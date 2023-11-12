function Y = inv_chol(L)
% Matrix Inversion using Cholesky Decomposition

N = size(L, 1) ;
Y = zeros(N, N) ;
% Work with the upper triangular matrix
R = L' ;
% Construct the auxillary diagonal matrix S = 1/rii
S = inv(diag(diag(R))) ;

for j=N:-1:1
    for i=j:-1:1
        Y(i,j) = S(i,j) - R(i,i+1:end)*Y(i+1:end,j) ;
        Y(i,j) = Y(i,j)/R(i,i) ;
        % Write out the symmetric element
        Y(j,i) = conj(Y(i,j)) ;
    end
end
end