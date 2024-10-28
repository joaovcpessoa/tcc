M = 8;
K = 4;

% Matriz de canal H (M x K)
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);

% Geração de bits aleatórios b para modulação 4-QAM
N_bits = 2;
bits = randi([0 1], K*N_bits, 1);

% Separação dos bits em grupos de 2 para mapeamento
bit_pairs = bits(:, 1:2:end) * 2 + bits(:, 2:2:end);

% Mapeamento 4-QAM usando a função qammod
M_qam = 4; % Ordem da modulação (4-QAM)
symbols = qammod(bit_pairs,M_qam);

% Transmissão através do canal H
% Recebido: y = H * s
received_signal = H * symbols;

% AWGN com variância 10^-3
noise_variance = 10^-3;
noise = sqrt(noise_variance/2) * (randn(size(received_signal)) + 1i * randn(size(received_signal)));
received_signal_with_noise = received_signal + noise;

disp('Sinal recebido através do canal H (com ruído):');
disp(received_signal_with_noise);