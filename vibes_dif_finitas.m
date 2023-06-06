% Parâmetros do sistema
m = 1;     % massa (kg)
c = 5;     % constante do amortecedor (Ns/m)
k = 1000;  % constante da mola (N/m)

% Condições iniciais
x0 = 0.7;  % posição inicial (m)
v0 = 30;   % velocidade inicial (m/s)

% Função da força externa
F = @(t) -100 *cos(50*t);  % Exemplo: F(t) = -100*cos(50t)

% Tempo de simulação
t_start = 0;
t_end = 5;
dt = 0.02;
t = t_start:dt:t_end;


%Iterações a partir das condições iniciais
x(1) = x0;
x(2) = x(1) + v0*dt;


%Resto das iterações
for i = 3 : (length(t))
    x(i) = x(i-1)*(2-dt*c/m - (dt)^2*k/m)+x(i-2)*(dt*c/m-1)-(dt)^2/m*F(t_start+i*dt);
end


%TESTE:
%
%teste_funct = m*diff(diff(x))/dt^2 + c*diff(x)(1:length(x)-2)/dt + k*(x)(1:length(x)-2);
%plot(t(1:length(t)-2),abs(teste_funct-F(t)(1:length(t)-2)),'bo-');
%
%O plot acima mostra a diferença em módulo entre " mx'' + cx' + kx " e " F(t) "
%
%O código abaixo mostra os valores de x0 e v0 esperados assim como seus valores
%obtidos pela x(t) calculada.
%
%x0
%v0
%disp(["x(0) = "  num2str(x(1))])
%disp(["x'(0) = " , num2str((diff(x)/dt)(1))]);

%Gráfico
length(t);
length(x);
figure;
plot(t,x,'bo-');
xlabel('Tempo (s)');
ylabel('Posição (m)');
title('Resposta Total do sistema (Diferenças finitas)');
