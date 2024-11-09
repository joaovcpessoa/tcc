% This Matlab script can be used to reproduce ... in the monograph:
% Joao Vitor, Correia Pessoa (2024), 
% "Análise de Impacto do Front-end Analógico em SistemaS MASSIVE MIMo",
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.
%
% Creation Date: 2024-11-03
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

% ------------------------------------------------------------------------------ %
% MEMORY ALLOCATION
% ------------------------------------------------------------------------------ %
s = zeros(N_BLK, K);

% ------------------------------------------------------------------------------ %
% MODULATION AND PRE-PROCESSING
% ------------------------------------------------------------------------------ %

bit_array = randi([0,1], B*N_BLK, K);

for users_idx = 1:K
    s(:, users_idx) = qammod(bit_array(:,users_idx), M_QAM, 'InputType', 'bit');
end
     
P_ZF = conj(H) / (H.' * conj(H));
P_MF = conj(H);

x_ZF = P_ZF * s.';
x_MF = P_MF * s.';

y_ZF = H.'* x_ZF;
y_MF = H.'* x_MF;

% ------------------------------------------------------------------------------ %
% APLICAÇÃO DO CLIPPING NOS SINAIS TRANSMITIDOS E RECEBIDOS
% ------------------------------------------------------------------------------ %

noise_variance = 10^-3;
P_noise = noise_variance;

SNR_dB = 20;
SNR = 10^(SNR_dB/10);
P_signal = SNR * P_noise;

A0 = sqrt(2 * P_signal);
clip = @(x, A0) min(abs(x), A0) .* exp(1j * angle(x));


x_ZF_clipped = clip(x_ZF, A0);
x_MF_clipped = clip(x_MF, A0);

y_ZF_clipped = H.' * x_ZF_clipped;
y_MF_clipped = H.' * x_MF_clipped;

% ------------------------------------------------------------------------------ %
% CÁLCULO DE BER
% ------------------------------------------------------------------------------ %
% Demodulação e comparação de bits, similar ao que é feito sem clipping

% ZF com clipping
bit_received_ZF_clipped = zeros(B*N_BLK, K);
for users_idx = 1:K
    s_received_ZF_clipped = y_ZF_clipped(users_idx, :).';
    bit_received_ZF_clipped(:, users_idx) = qamdemod(s_received_ZF_clipped, M_QAM, 'OutputType', 'bit');
end

% MF com clipping
bit_received_MF_clipped = zeros(B*N_BLK, K);
for users_idx = 1:K
    s_received_MF_clipped = y_MF_clipped(users_idx, :).';
    bit_received_MF_clipped(:, users_idx) = qamdemod(s_received_MF_clipped, M_QAM, 'OutputType', 'bit');
end

% BER para ZF com clipping
bit_errors_ZF_clipped = sum(bit_array ~= bit_received_ZF_clipped, 'all');
BER_ZF_clipped = bit_errors_ZF_clipped / numel(bit_array);

% BER para MF com clipping
bit_errors_MF_clipped = sum(bit_array ~= bit_received_MF_clipped, 'all');
BER_MF_clipped = bit_errors_MF_clipped / numel(bit_array);

% Resultados com clipping
disp('ZF com Clipping:')
fprintf('Number of bit errors (Clipped): %d, BER: %.5f\n', bit_errors_ZF_clipped, BER_ZF_clipped);

disp('')
disp('--------------------------------------------')
disp('')

disp('MF com Clipping:')
fprintf('Number of bit errors (Clipped): %d, BER: %.5f\n', bit_errors_MF_clipped, BER_MF_clipped);
