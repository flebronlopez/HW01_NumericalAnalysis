%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universidad de Puerto Rico en Rio Piedras
% Francheska Lebron Lopez
%           METODO NEWTON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_k, fx_k, iter, n_f, n_df, hist] = newton_method2(f, f_der, x_0, tol, maxit)

x_k = x_0; 

n_f  = 0; % cuenta las veces que se evalua la funcion f
n_df = 0; % cuenta las veces que se evalua la derivada de la funcion f

hist = x_k;  % % guarda xi's

% evalua f
fx_k = f(x_k); 
n_f = n_f + 1;

for iter = 1:maxit

    dfx_k = f_der(x_k); % evalua la derivada de la funcion f
    n_df = n_df + 1;

    if abs(dfx_k) < eps % verifica que la derivada no sea muy cercana a cero
        warning('Derivada ~0. Metodo detenido.');
        break
    end

    x_k_1 = x_k - fx_k/dfx_k; % evalua el termino x_k+1

    % guardar iteracion
    hist(end+1,1) = x_k_1;

    % evaluar f en el nuevo punto
    fx_k_1 = f(x_k_1);
    n_f = n_f + 1;

    % criterios de parada
    if (abs(fx_k_1) < tol) || (abs(x_k_1 - x_k) < tol)
        x_k = x_k_1;
        fx_k = fx_k_1;
        break
    end

    x_k  = x_k_1; % actualiza  
    fx_k = fx_k_1;

end

end
