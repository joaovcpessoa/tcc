% This Matlab script can be used to reproduce ... in the monograph:
% Joao Vitor, Correia Pessoa (2024), 
% "Análise de Impacto do Front-end Analógico em SistemaS MASSIVE MIMo",
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.
%
% Creation Date: 2024-10-27
% Last Updated:  2024-11-03
% ------------------------------------------------------------------------------ %
% GENERAL PARAMETERS
% ------------------------------------------------------------------------------ %

M = 100;
K = 4;
N_BLK = 1;
B = 8;
M_QAM = 2^B;
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);

% Distorção AM-AM e AM-PM
chi_A = 1;
kappa_A = 0.25;
chi_phi = 0.26;
kappa_phi = 0.25;

% ------------------------------------------------------------------------------ %
% MEMORY ALLOCATION
% ------------------------------------------------------------------------------ %
s = zeros(N_BLK, K);

% ------------------------------------------------------------------------------ %
% MODULATION AND PRE-PROCESSING
% ------------------------------------------------------------------------------ %

bit_array = randi([0, 1], B*N_BLK, K);

for users_idx = 1:K
    s(:, users_idx) = qammod(bit_array(:, users_idx), M_QAM, 'InputType', 'bit');
end

P_ZF = conj(H) / (H.' * conj(H));
P_MF = conj(H) ./ (vecnorm(H).^2);

x_ZF = P_ZF * s.';
x_MF = P_MF * s.';

x_ZF_distorted = (chi_A * x_ZF) ./ (1 + kappa_A * abs(x_ZF).^2);
x_MF_distorted = (chi_A * x_MF) ./ (1 + kappa_A * abs(x_MF).^2);

y_ZF = H.' * x_ZF_distorted;
y_MF = H.' * x_MF_distorted;

% ------------------------------------------------------------------------------ %
% DEMODULATION AND BIT COMPARISON
% ------------------------------------------------------------------------------ %

% ZF
bit_received_ZF = zeros(B*N_BLK, K);
for users_idx = 1:K
    s_received_ZF = y_ZF(users_idx, :).';
    bit_received_ZF(:, users_idx) = qamdemod(s_received_ZF, M_QAM, 'OutputType', 'bit');
end

% MF
bit_received_MF = zeros(B*N_BLK, K); % Reset for MF
for users_idx = 1:K
    s_received_MF = y_MF(users_idx, :).';
    bit_received_MF(:, users_idx) = qamdemod(s_received_MF, M_QAM, 'OutputType', 'bit');
end

% ------------------------------------------------------------------------------ %
% BIT ERROR RATE (BER) CALCULATION
% ------------------------------------------------------------------------------ %

noise_variance = 10^-3;

noise_ZF = sqrt(noise_variance/2) * (randn(size(s_received_ZF)) + 1i * randn(size(s_received_ZF)));
noise_MF = sqrt(noise_variance/2) * (randn(size(s_received_MF)) + 1i * randn(size(s_received_MF)));

received_signal_with_noise_ZF = s_received_ZF + noise_ZF;
received_signal_with_noise_MF = s_received_MF + noise_MF;

% BER for ZF
bit_errors_ZF = sum(bit_array ~= bit_received_ZF, 'all');
BER_ZF = bit_errors_ZF / numel(bit_array);

% BER for MF
bit_errors_MF = sum(bit_array ~= bit_received_MF, 'all');
BER_MF = bit_errors_MF / numel(bit_array);

% Results
disp('Análise com Pré-codificador Zero Forcing:');
disp('--')
fprintf('Sinal recebido através do canal H: %.4f + %.4fi\n', real(s_received_ZF), imag(s_received_ZF));
fprintf('Sinal recebido através do canal H (com ruído): %.4f + %.4fi\n', real(received_signal_with_noise_ZF), imag(received_signal_with_noise_ZF));
fprintf('Número de bits com erro: %d\n', bit_errors_ZF);
fprintf('Taxa de erro de bit (BER): %.5f\n', BER_ZF);
fprintf('\n')

disp('Análise com Pré-codificador Matched Filter:')
disp('--')
fprintf('Sinal recebido através do canal H: %.4f + %.4fi\n', real(s_received_MF), imag(s_received_MF));
fprintf('Sinal recebido através do canal H (com ruído): %.4f + %.4fi\n', real(received_signal_with_noise_MF), imag(received_signal_with_noise_MF));
fprintf('Número de bits com erro: %d\n', bit_errors_MF);
fprintf('Taxa de erro de bit (BER): %.5f\n', BER_MF);