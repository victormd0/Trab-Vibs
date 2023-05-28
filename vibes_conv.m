% Par�metros do sistema
m = 1;     % massa (kg)
c = 5;     % constante do amortecedor (Ns/m)
k = 1000;  % constante da mola (N/m)

% Condi��es iniciais
x0 = 0.7;  % posi��o inicial (m)
v0 = 30;   % velocidade inicial (m/s)

% Fun��o da for�a externa
F = @(t) -100 * (50 * t);  % Exemplo: F(t) = -100*(50t)

% Tempo de simula��o
t_start = 0;
t_end = 5;
dt = 0.02;
t = t_start:dt:t_end;

% Frequ�ncia natural do sistema
omega_n = sqrt(k / m);

% Fator de amortecimento
xi = c / (2 * sqrt(m * k));

% Resposta homog�nea
omega_d = omega_n * sqrt(1 - xi^2);
c1 = x0;
c2 = (v0 + xi * x0 * omega_n) / omega_d;
x_hom = @(t) exp(-xi * omega_n * t) .* (c1 * cos(omega_d * t) + c2 * sin(omega_d * t));

% C�lculo da resposta total pela integral de convolu��o
response = zeros(size(t));
for i = 1:length(t)
    integral = 0;
    for j = 1:i
        integral = integral + x_hom(t(i)-t(j)) * F(t(i)-t(j)) * dt;
    end
    response(i) = integral;
end

% Plot da resposta total
figure;
plot(t, response);
xlabel('Tempo (s)');
ylabel('Posi��o (m)');
title('Resposta Total do sistema MMA (Integral de Convolu��o)');