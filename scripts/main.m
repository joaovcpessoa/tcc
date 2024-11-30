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

savefig = 1;
addpath('./functions/');
root_save = 'C:\Users\joaov_zm1q2wh\OneDrive\Code\github\tcc\images\';

% ####################################################################### %
% PARÂMETROS PRINCIPAIS
% ####################################################################### %

N_BLK = 10000;

M = 100;
K = 4;
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);

B = 4;
M_QAM = 2^B;

SNR = -10:20;
N_SNR = length(SNR);
snr = 10.^(SNR/10);

N_PRE = {'MF', 'ZF', 'MMSE'};
N_AMP = {'IDEAL', 'CLIPPING', 'TWT', 'SOLID_STATE'};

BER = zeros(K, N_SNR, N_PRE, N_AMP);

A0 = 1;  % Fator de recorte
rho = 1; % Fator de suavidade da função não linear

% ####################################################################### %
% MODELO DO CANAL E DEFINIÇÃO DOS PRECODIFICADORES
% ####################################################################### %

bit_array = randi([0, 1], B * N_BLK, K);

% Debbugar para tirar prova real de que s=s1
% s1 = zeros(N_BLK, K);
% for users_idx = 1:K
%     s1(:, users_idx) = qammod(bit_array(:, users_idx), M_QAM, 'InputType', 'bit'); 
% end
s = qammod(bit_array, M_QAM, 'InputType', 'bit');

Ps = vecnorm(s).^2 / N_BLK;

v = sqrt(0.5) * (randn(K, N_BLK) + 1i*randn(K, N_BLK));
Pv = vecnorm(v,2,2).^2/N_BLK;
v_normalized = v./sqrt(Pv);

% ####################################################################### %
% TRANSMISSÃO E RECEPÇÃO
% ####################################################################### %

% Loop principal para SNR e amplificadores
for snr_idx = 1:N_SNR
    for pre_idx = 1:N_PRE
        for amp_idx = 1:N_AMP
           y = precoder(x, type);

            bit_received = zeros(B*N_BLK, K);

            for users_idx = 1:K
                s_received = y(users_idx, :).';
                Ps_received = norm(s_received)^2/N_BLK;
                bit_received(:, users_idx) = qamdemod(sqrt(Ps(users_idx)/Ps_received) * s_received, M_QAM, 'OutputType', 'bit');
                [~, BER(users_idx, snr_idx, pre_idx, amp_idx)] = biterr(bit_received, bit_array);
            end
        end
    end
end

% ####################################################################### %
% PLOTAGEM
% ####################################################################### %

figure;
set(gcf,'position',[0 0 800 600]);

for pre_idx = 1:N_PRE
    for amp_idx = 1:N_AMP
        semilogy(SNR, BER_ZF(:,pre_idx),  'o-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(pre_idx,:)); 
        hold on;
        % semilogy(SNR, BER_MF(:,pre_idx),  's-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(pre_idx+1,:));
        % semilogy(SNR, BER_MMSE(:,pre_idx),'v-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(pre_idx+2,:));
end

xlabel('SNR (dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);
legend({'ZF_{IDEAL}', 'ZF_{CLIPPING_A0}', 'ZF_{CLIPPING_A1}', 'ZF_{CLIPPING_A2}', 'ZF_{CLIPPING_A3}', 'ZF_{CLIPPING_A4}', 'ZF_{CLIPPING_A5}'}, 'Location', 'northeast', 'FontSize', fontsize);
legend box off;

set(gca, 'FontName', fontname, 'FontSize', fontsize);

if savefig == 1
    saveas(gcf,[root_save 'ber'],'fig');
    saveas(gcf,[root_save 'ber'],'png');
    saveas(gcf,[root_save 'ber'],'epsc2');
end