% Par�metros do sistema
m = 1;     % massa (kg)
c = 5;     % constante do amortecedor (Ns/m)
k = 1000;  % constante da mola (N/m)


% Condi��es iniciais
x0 = 0.7;  % posi��o inicial 0.7(m)
v0 = 30;   % velocidade inicial 30(m/s)

% Fun��o da For�a Externa
F = @(t) -100 * cos(50 * t);

% Tempo de simula��o
t_start = 0;
t_end = 5;
dt = 0.02;
t = t_start:dt:t_end;

% Frequ�ncia natural do sistema
omega_n = sqrt(k / m);

% Fator de amortecimento
xi = c / (2 * sqrt(m * k));


omega_d = omega_n * sqrt(1 - xi^2);

Green = @(t) exp(-xi*omega_n*t)/(m*omega_d).*sin(omega_d*t);

% C�lculo da resposta parcial pela integral de convolu��o

xp = dt*conv(Green(t),F(t));
xp = xp(1:length(t));


% C�lculo de c1 e c2
%diff: Fun��o para calcular diferen�as forward e, com isso, calcular a derivada

c1 = x0 - xp(1);
derxp = diff(xp)/dt;
c2 = (v0 - derxp(1)+xi*omega_n*c1)/omega_d;
x_hom = @(t) exp(-xi * omega_n * t) .* (c1 * cos(omega_d * t) + c2 * sin(omega_d * t));
x = x_hom(t) + xp;

%TESTE:
%
%teste_funct = m*diff(diff(x))/dt^2 + c*diff(x)(1:length(x)-2)/dt + k*(x)(1:length(x)-2);
%plot(t(1:length(t)-2),abs(teste_funct-F(t)(1:length(t)-2)),'bo-');
%
%O plot acima mostra a diferen�a em m�dulo entre " mx'' + cx' + kx " e " F(t) "
%Observa-se que o plot diminui em valor conforme dt diminui
%Se dt = 0.02, o valor m�ximo do plot � por volta de 530
%Se dt = 0.0001, o valor m�ximo do plot � por volta de 2.8
%
%O c�digo abaixo mostra os valores de x0 e v0 esperados assim como seus valores
%obtidos pela x(t) calculada.
%
%x0
%v0
%disp(["x(0) = "  num2str(x(1))])
%disp(["x'(0) = " , num2str((diff(x)/dt)(1))]);
%
%Novamente, observa-se uma converg�ncia de x(0) e x'(0) para x0 e v0 conforme dt diminui
%


% Plot da resposta total
figure;
plot(t,x,'bx-');
xlabel('Tempo (s)');
ylabel('Posi��o (m)');
title('Resposta Total do sistema (Integral de Convolu��o)');

data = [x' t'];
csvwrite('C:\Users\Boiling\Desktop\vis comp\DADOS_DO_MATLAB.csv', data);