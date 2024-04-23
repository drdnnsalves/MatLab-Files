%Input variables:
%  A: Square matrix;
%  b: collumn vector;
%  x_0: Initial aproximation of the system;
%  E: Tolerency;
%  M: Max number of iterations;
%  N: Norm (1, 2, inf)
%
%Output variables
%  x_k: Solution found;
%  n: Norm of the last 2 iterations;
%  k: Number of iterations performed;
%  r: Norm of the residual

function [x_k, n, k, r] = gauss_seidel_2(A, b, x_0, E, M, N)

    % Separamos a matriz A em superior (U), diagonal (D) e inferior (L)
    L = tril(A, -1);
    D = diag(diag(A));
    U = triu(A, 1);
     
    % Inicializando o vetor x_k e o contador k
    x_k = x_0; 
    k = 0;
    
    % Tomando o valor inicial para n
    n = norm(b - A*x_0, N);

    while  k < M && n > E
        % Guardando o valor da iteração anterior
        x_ant = x_k;
        
        % Iterando o valor de x_k usando a função solve_system_inf do
        % arquivo solve_system_inf.m
        x_k = solve_system_inf(L,D,U,x_ant,b);

        n = norm(x_k - x_ant, N);
        k = k + 1;
      
    end
    % Norma do resíduo feita após as iterações
    r = norm(b - A*x_k, N);

end
