%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Universidad de Puerto Rico en Rio Piedras
% Francheska Lebron Lopez
%           METODO BISECCION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xk, fxk, iter, n_f, hist] = bisection_method2(f, a, b, tol, maxit)

n_f = 0; % cuenta las veces que se evalua la funcion f

fa = f(a); n_f = n_f + 1; % evalua la funcion en el extremo a
fb = f(b); n_f = n_f + 1; % evalua la funcion en el extremo b

if fa * fb > 0 % verifica si f(a)f(b) es negativo
    error('El intervalo no satisface f(a)f(b) < 0.');
end

hist = []; % guarda aproximaciones

for iter = 1:maxit

    c = (a + b)/2; % calcula el valor medio del intervalo
    fc = f(c); n_f = n_f + 1; % evalua la funcion en c

    hist(end+1,1) = c;

    % criterios de parada
    if (abs(fc) < tol) || (abs(b - a)/2 < tol)
        break
    end

    if fa * fc < 0 % si f(a)f(c) son  negativos elige el valor c para el nuevo intervalo
        b = c; fb = fc;
    else
        a = c; fa = fc;
    end

end

xk  = c; % actualiza
fxk = fc;

end
