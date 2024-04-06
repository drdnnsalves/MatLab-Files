%  Jacobi Method for resolution of linear systems.
%
% Imput Variables:
%   A: Square matrix;
%   b: Collumn vector;
%   x_0: Initial aproximation to the system's solution;
%   E: Tolerance;
%   M: Maximum number of iterations;
%   N: Norm type (Possible values: 1, 2 or inf).
% 
% Output:
%   x_k: Solution; 
%   n: Norm of the difference of the last two iterations;
%   k: Number of iterations performed;
%   r: Residual's norm.

function [x_k, n, k, r] = jacobi_method_1(A, b, x_0, E, M, N)

    % We separate the matrix A in Upper (U), Diagonal (D) and Lower (L) 
    L = tril(A, -1);
    D = diag(diag(A));
    U = triu(A, 1);
     
    % Initializing the vector x_k and the counter k
    x_k = x_0; 
    k = 0;
    
    % Taking the initial value for n
    n = norm(b - A*x_0, N);

    while  k < M && n > E
        % Keeping the value of the last iteration
        x_old = x_k;
        
        x_k = (b - (L + U) * x_old) ./ diag(D);
        
        n = norm(x_k - x_old, N);
        k = k + 1;
      
    end
    % Residual norm after the iterations
    r = norm(b - A*x_k, N);

end
