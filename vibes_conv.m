% Parâmetros do sistema
m = 1;     % massa (kg)
c = 5;     % constante do amortecedor (Ns/m)
k = 1000;  % constante da mola (N/m)

% Condições iniciais
x0 = 0.7;  % posição inicial (m)
v0 = 30;   % velocidade inicial (m/s)

% Função da força externa
F = @(t) -100 * (50 * t);  % Exemplo: F(t) = -100*(50t)

% Tempo de simulação
t_start = 0;
t_end = 5;
dt = 0.02;
t = t_start:dt:t_end;

% Frequência natural do sistema
omega_n = sqrt(k / m);

% Fator de amortecimento
xi = c / (2 * m);

% Resposta homogênea
omega_d = omega_n * sqrt(1 - xi^2);
c1 = x0;
c2 = (v0 + xi * x0 * omega_n) / omega_d;
x_hom = @(t) exp(-xi * omega_n * t) .* (c1 * cos(omega_d * t) + c2 * sin(omega_d * t));

% Resposta permanente
F0 = -100 * (50 / sqrt((k - m * omega_n^2)^2 + (c * omega_n)^2));
x_per = @(t) F0 * cos(omega_n * t - atan(c * omega_n / (k - m * omega_n^2)));

% Cálculo da resposta total pela integral de convolução usando a função conv
response = conv(x_hom(t), F(t), 'same') * dt + x_per(t);

% Plot da resposta total
figure;
plot(t, response);
xlabel('Tempo (s)');
ylabel('Posição (m)');
title('Resposta Total do sistema MMA (Integral de Convolução)');