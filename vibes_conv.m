% Parâmetros do sistema
m = 1;     % massa (kg)
k = 1000;  % constante da mola (N/m)
c = 50;    % constante do amortecedor (Ns/m)

% Força harmônica
F = @(t) -100*(50*t);  % F(t) = -100*(50t) (N)

% Condições iniciais
x0 = 0.7;    % posição inicial (m)
v0 = 30;    % velocidade inicial (m/s)

% Tempo de simulação
t_start = 0;
t_end = 5;
dt = 0.02;
t = t_start:dt:t_end;
n = length(t);

% Coeficientes da resposta homogênea e permanente
omega = sqrt(k/m - (c/(2*m))^2);
alpha = c/(2*m);

if c < 2*m*omega
    C1 = (v0 + alpha*x0)/omega;
    C2 = x0;
else
    alpha_1 = (-c + sqrt(c^2 - 4*m*k))/(2*m);
    alpha_2 = (-c - sqrt(c^2 - 4*m*k))/(2*m);
    C1 = (v0 - alpha_2*x0)/(alpha_1 - alpha_2);
    C2 = x0 - C1;
end

% Cálculo da resposta total pela integral de convolução
response = zeros(1, n);

for i = 1:n
    x_hom = @(tau) (C1*cos(omega*(t(i)-tau)) + C2*sin(omega*(t(i)-tau))) .* exp(-alpha*(t(i)-tau));
    x_per = @(tau) F(tau) / k;
    response(i) = integral(@(tau) x_hom(tau) .* x_per(t(i)-tau), 0, t(i));
end

% Adicionando a resposta homogênea
response = response + x_hom(t);

% Plot da resposta
figure;
plot(t, response);
xlabel('Tempo (s)');
ylabel('Resposta (m)');
title('Resposta do sistema massa-mola-amortecedor pela integral de convolução');