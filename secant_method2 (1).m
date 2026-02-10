%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universidad de Puerto Rico en Rio Piedras
% Francheska Lebron Lopez
%           METODO SECANTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xk, fxk, iter, n_f, hist] = secant_method2(f, x0, x1, tol, maxit)

n_f = 0; % cuenta las veces que se evalua la funcion f

fx0 = f(x0); n_f = n_f + 1; % evalua f(x_0)
fx1 = f(x1); n_f = n_f + 1; % evalua f(x_1)

hist = [x0; x1]; % guarda xi's

% inicializa por seguridad
x2 = x1;
fx2 = fx1;

for iter = 1:maxit

    denominador = (fx1 - fx0); % computa el denominador

    % verifica qe el denominador no sea un numero muy cercano a cero
    if abs(denominador) < 1e-14 * max(1, abs(fx1) + abs(fx0))
        warning('Metodo detenido: f(x1) - f(x0) ~ 0.');
        x2 = x1;
        fx2 = fx1;
        break
    end

    x2 = x1 - fx1 * (x1 - x0) / denominador; % computa x_2

    fx2 = f(x2); %evalua la funcion en el nuevo punto
    n_f = n_f + 1;

    hist(end+1,1) = x2;

    % criterios de parada
    if (abs(fx2) < tol) || (abs(x2 - x1) < tol)
        break
    end
    % actualiza para nuevo intervalo
    x0 = x1; f 
    x0 = fx1; 
    x1 = x2; fx1 = fx2;

end

xk  = x2;
fxk = fx2;

end
