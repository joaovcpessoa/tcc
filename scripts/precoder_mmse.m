clear;
close all;
clc;

% ####################################################################### %
%% PARÂMETROS PRINCIPAIS
% ####################################################################### %

addpath('./functions/');

N_BLK = 10000;

M = 100;
K = 4;
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);

B = 4;
M_QAM = 2^B;

SNR = -10:50;
N_SNR = length(SNR);
snr = 10.^(SNR/10);

N_A0 = 5;
N_AMP = 4;

A0 = [0.5, 1.0, 1.5, 2.0, 2.5];

precoder_type = 'MMSE';
amplifiers_type = {'IDEAL', 'CLIP', 'TWT', 'SS'};

BER = zeros(K, N_SNR, N_AMP, N_A0);

% ####################################################################### %
%% PARAMETROS DO SINAL E DO RUÍDO 
% ####################################################################### %

bit_array = randi([0, 1], B * N_BLK, K);
s = qammod(bit_array, M_QAM, 'InputType', 'bit');
Ps = vecnorm(s).^2 / N_BLK;

precoder = compute_precoder(precoder_type, H, N_SNR, snr);
x_normalized = normalize_precoded_signal(precoder, precoder_type, M, s, N_SNR);

v = sqrt(0.5) * (randn(K, N_BLK) + 1i*randn(K, N_BLK));
Pv = vecnorm(v,2,2).^2/N_BLK;
v_normalized = v./sqrt(Pv);

% ####################################################################### %
%% TRANSMISSÃO E RECEPÇÃO
% ####################################################################### %

for snr_idx = 1:N_SNR
    for a_idx = 1:N_A0
        for amp_idx = 1:N_AMP
            a0 = A0(a_idx);
            current_amp_type = amplifiers_type{amp_idx};

            y = H.' * amplifier(sqrt(snr(snr_idx)) * x_normalized(:, :, snr_idx), current_amp_type, a0) + v_normalized;
            bit_received = zeros(B*N_BLK, K);

            for users_idx = 1:K
                s_received = y(users_idx, :).';
                Ps_received = norm(s_received)^2/N_BLK;
                bit_received(:, users_idx) = qamdemod(sqrt(Ps(users_idx)/Ps_received) * s_received, M_QAM, 'OutputType', 'bit');
                [~, BER(users_idx, snr_idx, amp_idx, a_idx)] = biterr(bit_received(:, users_idx), bit_array(:, users_idx));
            end
        end
    end
end

save('ber_mmse.mat','BER','y','SNR', 'N_AMP', 'N_A0');