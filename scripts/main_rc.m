% ####################################################################### %
% LIMPEZA
% ####################################################################### %

clear;
close all;
clc;

% ####################################################################### %
% PARÂMETROS DE PLOTAGEM
% ####################################################################### %

linewidth  = 2;
fontname   = 'Times New Roman';
fontsize   = 20;
markersize = 10;

colors = [0.0000 0.0000 0.0000;
          0.0000 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.9290 0.6940 0.1250;
          0.4940 0.1840 0.5560;
          0.4660 0.6740 0.1880;
          0.3010 0.7450 0.9330;
          0.6350 0.0780 0.1840];

savefig   = 1;
root_save = 'C:\Users\joaov_zm1q2wh\OneDrive\Code\github\tcc\images\';

% ####################################################################### %
% PARÂMETROS PRINCIPAIS
% ####################################################################### %

N_BLK = 10;

M = 100;
K = 4;
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);

B = 4;
M_QAM = 2^B;

SNR = -10:20;
N_SNR = length(SNR);
snr = 10.^(SNR/10);

N_PRE = 4;
BER_ZF = zeros(1, N_SNR, N_PRE);
BER_MF = zeros(1, N_SNR, N_PRE);
BER_MMSE = zeros(1, N_SNR, N_PRE);

% ####################################################################### %
% PARÂMETROS DE MODELAGEM DAS NÃO LINEARIDADES
% ####################################################################### %

A0 = 1;  % Fator de recorte
rho = 1; % Fator de suavidade da função não linear

% Amplificador com corte ideal
ideal_clipping_amp = @(x, A0) min(abs(x), A0) .* exp(1j * angle(x));

% Amplificador TWT
chi_A = 1;
kappa_A = 0.25;
chi_phi = 0.26;
kappa_phi = 0.25;
twt_amp1 = @(x) abs(x) .* (1 - kappa_A .* abs(x).^2) .* exp(1j * (angle(x) + chi_A .* abs(x).^2 + kappa_A .* abs(x).^4));       % AM-AM
twt_amp2 = @(x) abs(x) .* (1 - kappa_phi .* abs(x).^2) .* exp(1j * (angle(x) + chi_phi .* abs(x).^2 + kappa_phi .* abs(x).^4)); % AM-PM

% Amplificador de estado sólido
solid_state_amp = @(x, A0, rho) abs(x) ./ (1 + (abs(x) / A0).^(2 * rho)).^(1 / (2 * rho)) .* exp(1j * angle(x));

% Polinomial
AIIP3 = (4/3 * kappa_A)^(1/2);
AIIP5 = (8/5 * kappa_A^2)^(1/4);
% AIIP3 = (8/3)^(1/2) * A0;
% AIIP5 = (64/3)^(1/4) * A0;
b3 = 4/(3*AIIP3^2);
b5 = 8/(5*AIIP5^4);
polinomial = @(x) abs(x) .* (1 + b3 * abs(x).^2 + b5 * abs(x).^4) .* exp(1j * angle(x));

% Amplificadores: no_amp, ideal_clipping_amp, twt_amp, solid_state_amp
amplifiers = {@(x) x, @(x) ideal_clipping_amp(x, A0), twt_amp1, @(x) solid_state_amp(x, A0, rho), polinomial};

% ####################################################################### %
% MODELO DO CANAL E DEFINIÇÃO DOS PRECODIFICADORES
% ####################################################################### %

s = zeros(N_BLK, K);
bit_array = randi([0, 1], B * N_BLK, K);
for users_idx = 1:K
    s(:, users_idx) = qammod(bit_array(:, users_idx), M_QAM, 'InputType', 'bit'); 
end

Ps = vecnorm(s).^2 / N_BLK;

% Pré-codificadores
precoder_ZF = conj(H) / (H.' * conj(H));
precoder_MF = conj(H) ./ (vecnorm(H).^2);
precoder_MMSE = zeros(M, K, N_SNR);
for snr_idx = 1:N_SNR
    precoder_MMSE(:,:,snr_idx) = conj(H) / (H.' * conj(H) + 1/snr(snr_idx) * eye(K));
end

% ####################################################################### %
% TRANSMISSÃO E RECEPÇÃO
% ####################################################################### %

v = sqrt(0.5) * (randn(K, N_BLK) + 1i*randn(K, N_BLK));
Pv = vecnorm(v,2,2).^2/N_BLK;
v_normalized = v./sqrt(Pv);

% Loop principal para SNR e amplificadores
for snr_idx = 1:N_SNR
    for pre_idx = 1:N_PRE
        % ZF
        y_ZF = H.' * amplifiers{pre_idx}(sqrt(snr(snr_idx)) * precoder_ZF * s.') + v_normalized; 
        y_MF = H.' * amplifiers{pre_idx}(sqrt(snr(snr_idx)) * precoder_MF * s.') + v_normalized;
        y_MMSE = H.' * amplifiers{pre_idx}(sqrt(snr(snr_idx)) * precoder_MMSE(:,:,snr_idx) * s.') + v_normalized; 
        
        % Demodulação e BER para ZF
        bit_received_ZF = zeros(B*N_BLK, K);
        for users_idx = 1:K
            s_received_ZF = y_ZF(users_idx, :).';
            Ps_received_ZF = norm(s_received_ZF)^2/N_BLK;
            bit_received_ZF(:, users_idx) = qamdemod(sqrt(Ps(users_idx)/Ps_received_ZF) * s_received_ZF, M_QAM, 'OutputType', 'bit');
        end
        [~, BER_ZF(snr_idx, pre_idx)] = biterr(bit_received_ZF, bit_array);

        % Demodulação e BER para MF
        bit_received_MF = zeros(B*N_BLK, K);
        for users_idx = 1:K
            s_received_MF = y_MF(users_idx, :).';
            Ps_received_MF = norm(s_received_MF)^2/N_BLK;
            bit_received_MF(:, users_idx) = qamdemod(sqrt(Ps(users_idx)/Ps_received_MF) * s_received_MF, M_QAM, 'OutputType', 'bit');
        end
        [~, BER_MF(snr_idx, pre_idx)] = biterr(bit_received_MF, bit_array);

        % Demodulação e BER para MMSE
        bit_received_MMSE = zeros(B*N_BLK, K);
        for users_idx = 1:K
            s_received_MMSE = y_MMSE(users_idx, :).';
            Ps_received_MMSE = norm(s_received_MMSE)^2/N_BLK;
            bit_received_MMSE(:, users_idx) = qamdemod(sqrt(Ps(users_idx)/Ps_received_MMSE) * s_received_MMSE, M_QAM, 'OutputType', 'bit');
        end
        [~, BER_MMSE(snr_idx, pre_idx)] = biterr(bit_received_MMSE, bit_array);
    end
end

% ####################################################################### %
% PLOTAGEM
% ####################################################################### %

figure;
set(gcf,'position',[0 0 800 600]);

for pre_idx = 1:N_PRE
    semilogy(SNR, BER_ZF(:,pre_idx),  'o-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(pre_idx,:)); 
    hold on;
    semilogy(SNR, BER_MF(:,pre_idx),  's-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(pre_idx+1,:));
    semilogy(SNR, BER_MMSE(:,pre_idx),'v-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(pre_idx+2,:));
end

xlabel('SNR (dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);
legend({'ZF-1', 'ZF-2', 'ZF-3', 'ZF-4', 'MF-1', 'MF-2', 'MF-3', 'MF-4', 'MMSE-1', 'MMSE-2', 'MMSE-3', 'MMSE-4'}, 'Location', 'northeast', 'FontSize', fontsize);
legend box off

set(gca, 'FontName', fontname, 'FontSize', fontsize);