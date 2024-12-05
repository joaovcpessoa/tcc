clear;
close all;
clc;

% ####################################################################### %
%% PARÂMETROS PRINCIPAIS
% ####################################################################### %

addpath('./functions/');

N_BLK = 10000;    % Número de blocos
M = 100;          % Número de antenas de transmissão
K = 4;            % Número de usuários
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);  % Canal

B = 4;            % Número de bits por símbolo (para QAM)
M_QAM = 2^B;      % Ordem do QAM

SNR = -10:50;     % Intervalo de SNR (em dB)
N_SNR = length(SNR);
snr = 10.^(SNR/10);

N_A0 = 5;         % Número de valores de A0
N_AMP = 4;        % Número de tipos de amplificadores

A0 = [0.5, 1.0, 1.5, 2.0, 2.5];  % Valores de A0

precoder_type = 'MF';
amplifiers_type = {'IDEAL', 'CLIP', 'TWT', 'SS'};

BER = zeros(K, N_SNR, N_AMP, N_A0);  % Resultados BER
N_MC = 100;   % Número de simulações de Monte Carlo

% ####################################################################### %
%% PARÂMETROS DO SINAL E DO RUÍDO 
% ####################################################################### %

bit_array = randi([0, 1], B * N_BLK, K);  % Gerar os bits de entrada
s = qammod(bit_array, M_QAM, 'InputType', 'bit');  % Modulação QAM
Ps = vecnorm(s).^2 / N_BLK;  % Potência do sinal

precoder = compute_precoder(precoder_type, H, N_SNR, snr);
x_normalized = normalize_precoded_signal(precoder, precoder_type, M, s, N_SNR);

v = sqrt(0.5) * (randn(K, N_BLK) + 1i*randn(K, N_BLK));  % Ruído
Pv = vecnorm(v,2,2).^2 / N_BLK;  % Potência do ruído
v_normalized = v ./ sqrt(Pv);  % Normalização do ruído

% ####################################################################### %
%% SIMULAÇÃO DE MONTE CARLO
% ####################################################################### %

for mc_idx = 1:N_MC
    for snr_idx = 1:N_SNR
        for a_idx = 1:N_A0
            for amp_idx = 1:N_AMP
                a0 = A0(a_idx);
                current_amp_type = amplifiers_type{amp_idx};

                y = H.' * amplifier(sqrt(snr(snr_idx)) * x_normalized, current_amp_type, a0) + v_normalized;
                bit_received = zeros(B * N_BLK, K);

                for users_idx = 1:K
                    s_received = y(users_idx, :).';
                    Ps_received = norm(s_received)^2 / N_BLK;
                    bit_received(:, users_idx) = qamdemod(sqrt(Ps(users_idx) / Ps_received) * s_received, M_QAM, 'OutputType', 'bit');
                    
                    [~, bit_error] = biterr(bit_received(:, users_idx), bit_array(:, users_idx));
                    BER(users_idx, snr_idx, amp_idx, a_idx) = BER(users_idx, snr_idx, amp_idx, a_idx) + bit_error;
                end
            end
        end
    end
end

BER = BER / N_MC;

save('ber_mc_mf.mat', 'BER', 'y', 'SNR', 'N_AMP', 'N_A0');
