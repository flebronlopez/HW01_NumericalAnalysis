%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Universidad de Puerto Rico en Rio Piedras
% Francheska Lebron Lopez
% Num. Est. 842094050
% MAIN: Newton / Secante / Biseccion / Analisis de Convergencia
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%
% Funcion f(x) y f'(x)
%%%%%%%%%%%%%%%%%%%%%%
f  = @(x) x .* (x.^2 - 1) .* (x - 3) .* exp(-0.5 * (x - 1).^2);
fd = @(x) (3 + x - 13.*x.^2 + 2.*x.^3 + 4.*x.^4 - x.^5) .* exp(-0.5 * (x - 1).^2);

%%%%%%%%%%%%%%%%%%%%%%
% Parametros Generales
%%%%%%%%%%%%%%%%%%%%%%
rng(89);
tol   = 1e-15;
maxit = 10000;

ValorTeorico = [-1; 0; 1; 3];
n = length(ValorTeorico);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniciales Newton 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = -1 + 0.6*(2*rand-1);
x2 =  0 + 0.6*(2*rand-1);
x3 =  1 + 0.6*(2*rand-1);
x4 =  3 + 0.6*(2*rand-1);
x0 = [x1,x2,x3,x4];

%%%%%%%%%%%%%%%%%%%%%%
% Iniciales Secante 
%%%%%%%%%%%%%%%%%%%%%%
aS = [-1.5; -0.5;  0.5;  2.5];
bS = [-0.5;  0.5;  1.5;  3.5];

x0s = zeros(n,1);
x1s = zeros(n,1);
for k = 1:n
    x0s(k) = aS(k) + (bS(k)-aS(k))*rand;
    x1s(k) = aS(k) + (bS(k)-aS(k))*rand;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iniciales Biseccion (intervalos con cambio de signo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aB = zeros(n,1);
bB = zeros(n,1);

for k = 1:n
    ok = false;
    while ~ok
        if k == 1
            a = -1.5 + 0.4*(2*rand-1);
            b = -0.5 + 0.4*(2*rand-1);
        elseif k == 2
            a = -0.5 + 0.4*(2*rand-1);
            b =  0.5 + 0.4*(2*rand-1);
        elseif k == 3
            a =  0.5 + 0.4*(2*rand-1);
            b =  1.5 + 0.4*(2*rand-1);
        else
            a =  2.5 + 0.4*(2*rand-1);
            b =  3.5 + 0.4*(2*rand-1);
        end

        if a > b, tmp=a; a=b; b=tmp; end

        if f(a)*f(b) < 0
            ok = true;
            aB(k) = a;
            bB(k) = b;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CORRER METODOS + guardar historiales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aprox_N = zeros(n,1); fAprox_N = zeros(n,1); Iter_N = zeros(n,1);
EvalF_N = zeros(n,1); EvalFD_N = zeros(n,1);
Hist_N  = cell(n,1);

Aprox_S = zeros(n,1); fAprox_S = zeros(n,1); Iter_S = zeros(n,1);
EvalF_S = zeros(n,1);
Hist_S  = cell(n,1);

Aprox_B = zeros(n,1); fAprox_B = zeros(n,1); Iter_B = zeros(n,1);
EvalF_B = zeros(n,1);
Hist_B  = cell(n,1);

for k = 1:n
    [Aprox_N(k), fAprox_N(k), Iter_N(k), EvalF_N(k), EvalFD_N(k), Hist_N{k}] = ...
        newton_method2(f, fd, x0(k), tol, maxit);

    [Aprox_S(k), fAprox_S(k), Iter_S(k), EvalF_S(k), Hist_S{k}] = ...
        secant_method2(f, x0s(k), x1s(k), tol, maxit);

    [Aprox_B(k), fAprox_B(k), Iter_B(k), EvalF_B(k), Hist_B{k}] = ...
        bisection_method2(f, aB(k), bB(k), tol, maxit);
end

%%%%%%%%%%%%%%%%%%%%%%
% TABLAS (resumen)
%%%%%%%%%%%%%%%%%%%%%%
T_Newton = table(x0', ValorTeorico, Aprox_N, fAprox_N, abs(Aprox_N-ValorTeorico), Iter_N, EvalF_N, EvalFD_N, ...
    'VariableNames', {'ValInicial','ValTeo','Aprox','fAprox','Error','Iter','Eval_f','Eval_fder'});

T_Secante = table(x0s, x1s, ValorTeorico, Aprox_S, fAprox_S, abs(Aprox_S-ValorTeorico), Iter_S, EvalF_S, ...
    'VariableNames', {'x0','x1','ValTeo','Aprox','fAprox','Error','Iter','Eval_f'});

T_Biseccion = table(aB, bB, ValorTeorico, Aprox_B, fAprox_B, abs(Aprox_B-ValorTeorico), Iter_B, EvalF_B, ...
    'VariableNames', {'a','b','ValTeo','Aprox','fAprox','Error','Iter','Eval_f'});

disp('==================== METODO DE NEWTON ===================='), disp(T_Newton)
disp('==================== METODO DE SECANTE ===================='), disp(T_Secante)
disp('==================== METODO DE BISECCION ===================='), disp(T_Biseccion)

%%%%%%%%%%%%%%%%%%%%%%
% Analisis de Convergencia
% Metodos y raices
%%%%%%%%%%%%%%%%%%%%%%
for k = 1:n

    xi = ValorTeorico(k);

    errN = abs(Hist_N{k} - xi); errN(errN==0) = realmin;
    errS = abs(Hist_S{k} - xi); errS(errS==0) = realmin;
    errB = abs(Hist_B{k} - xi); errB(errB==0) = realmin;

    figure
    semilogy(0:length(errN)-1, errN, '-', 'LineWidth', 2); hold on
    semilogy(0:length(errS)-1, errS, '-', 'LineWidth', 2);
    semilogy(0:length(errB)-1, errB, '-', 'LineWidth', 2);

    grid on
    xlabel('Iteracion k')
    ylabel('|x_k - \xi|')
    title(['Comparacion de convergencia (raiz \xi = ', num2str(xi), ')'])
    legend('Newton','Secante','Biseccion','Location','southwest')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q_obs: Newton y Secante y Biseccion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qObs_N = NaN(n,1);
qObs_S = NaN(n,1);
qObs_B = NaN(n,1);

for k = 1:n

    xi = ValorTeorico(k);

    % Newton
    e = abs(Hist_N{k} - xi); e(e==0)=realmin;
    if length(e) >= 4
        qk = zeros(length(e)-2,1);
        for j = 2:length(e)-1
            qk(j-1) = log(e(j+1)/e(j)) / log(e(j)/e(j-1));
        end
        t = min(3,length(qk));
        qObs_N(k) = mean(qk(end-t+1:end));
    end

    % Secante
    e = abs(Hist_S{k} - xi); e(e==0)=realmin;
    if length(e) >= 4
        qk = zeros(length(e)-2,1);
        for j = 2:length(e)-1
            qk(j-1) = log(e(j+1)/e(j)) / log(e(j)/e(j-1));
        end
        t = min(3,length(qk));
        qObs_S(k) = mean(qk(end-t+1:end));
    end

    % Biseccion 
    e = abs(Hist_B{k} - xi); e(e==0)=realmin;
    if length(e) >= 4
        qk = zeros(length(e)-2,1);
        for j = 2:length(e)-1
            qk(j-1) = log(e(j+1)/e(j)) / log(e(j)/e(j-1));
        end
        t = min(3,length(qk));
        qObs_B(k) = mean(qk(end-t+1:end));
    end

end

T_q = table(ValorTeorico, qObs_N, qObs_S, qObs_B, ...
    'VariableNames', {'Raiz','qObs_Newton','qObs_Secante','qObs_Biseccion'});

disp('==================== ORDEN OBSERVADO q_obs ====================')
disp(T_q)
