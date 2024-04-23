%  Auxiliary function for Gauss Seidel method
%
%Input variables:
%  L: Lower matrix;
%  D: Diagonal matrix;
%  U: Upper matrix;
%  x: vector for Gauss Seidel iteration;
%  b: idem.
%
%Output variables
%  x_k: Solution of (L+D)*x_k = -U*x + b.
%
%

function [x_k] = solve_system_inf(L, D, U, x, b)

    % Inicializando constantes e nomeando k, T
    n = length(x);
    k = -U*x + b;
    T = L + D;
    x_k = zeros(n, 1);

    % Solução do sistema que queremos
    for i = 1:n
        x_k(i) = k(i);
        for j = 1:i-1
            x_k(i) = x_k(i) - T(i,j) * x_k(j);
        end
        x_k(i) = x_k(i) / T(i,i);
    end
    
end
