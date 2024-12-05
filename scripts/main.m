% ####################################################################### %
%% LIMPEZA
% ####################################################################### %

clear;
%close all;
clc;

% ####################################################################### %
%% PARÂMETROS DE PLOTAGEM
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

savefig = 0;
addpath('./functions/');
root_save = 'C:\Users\joaov_zm1q2wh\OneDrive\Code\github\tcc\images\';

% ####################################################################### %
%% PARÂMETROS PRINCIPAIS
% ####################################################################### %

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

% ####################################################################### %
%% ESCOLHA DO PRECODIFICADOR E MODELO DE NÃO LINEARIDADE
% ####################################################################### %

disp('Precoders:');
disp('1 - ZF');
disp('2 - MF');
disp('3 - MMSE');

choice1 = input('Choose precoder: ');

if choice1 == 1
    precoder_type = 'ZF';
elseif choice1 == 2
    precoder_type = 'MF';
elseif choice1 == 3
    precoder_type = 'MMSE';
else
    error('Invalid precoder type.');
end

disp('####################');

disp('Amplifiers:');
disp('1 - IDEAL');
disp('2 - CLIP');
disp('3 - TWT');
disp('4 - SS');

choice2 = input('Choose amplifiers: ');

if choice2 == 1
    amplifier_type = 'IDEAL';
elseif choice2 == 2
    amplifier_type = 'CLIP';
elseif choice2 == 3
    amplifier_type = 'TWT';
elseif choice2 == 4
    amplifier_type = 'SS';
else
    error('Invalid amplifier type.');
end

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

% Iterar sob o SNR
% Iterar sob o A
% Iterar sob o K

BER = zeros(K, N_SNR, N_AMP, N_A0);

for snr_idx = 1:N_SNR
    for a_idx = 1:N_A0
        switch upper(precoder_type)
            case {'ZF', 'MF'}
               y = H.' * amplifier(sqrt(snr(snr_idx)) * x_normalized, amplifier_type, a_idx) + v_normalized;
            case 'MMSE'
               y = H.' * amplifier(sqrt(snr(snr_idx)) * x_normalized(:, :, snr_idx), amplifier_type, a_idx) + v_normalized;
        end
        bit_received = zeros(B*N_BLK, K);

        for users_idx = 1:K
            s_received = y(users_idx, :).';
            Ps_received = norm(s_received)^2/N_BLK;
            bit_received(:, users_idx) = qamdemod(sqrt(Ps(users_idx)/Ps_received) * s_received, M_QAM, 'OutputType', 'bit');
            [~, BER(users_idx, snr_idx, a_idx)] = biterr(bit_received(:, users_idx), bit_array(:, users_idx));
        end
    end
end

save('ber_zf.mat','BER','y');

% ####################################################################### %
%% GRÁFICOS
% ####################################################################### %

figure;

set(gcf,'position',[0 0 800 600]);

for a_idx = 1:N_A0
    semilogy(SNR, mean(BER(:,:,a_idx),1), 'o-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(a_idx,:));
    hold on;
end

xlabel('SNR (dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);

legend({strcat(precoder_type, '-', amplifier_type, '_{A=1}'), ...
        strcat(precoder_type, '-', amplifier_type, '_{A=2}'), ...
        strcat(precoder_type, '-', amplifier_type, '_{A=3}'), ...
        strcat(precoder_type, '-', amplifier_type, '_{A=4}'), ...
        strcat(precoder_type, '-', amplifier_type, '_{A=5}')' ...
        }, 'Location', 'northeast', 'FontSize', fontsize);
legend box off;

set(gca, 'FontName', fontname, 'FontSize', fontsize);

graph_name = strcat('ber_', precoder_type, '_', amplifier_type);

if savefig == 1
    saveas(gcf,[root_save graph_name],'fig');
    saveas(gcf,[root_save graph_name],'png');
    saveas(gcf,[root_save graph_name],'epsc2');
end

% ####################################################################### %

%% DEBUG
% s_old = zeros(N_BLK, K);
% for users_idx = 1:K
%     s_old(:, users_idx) = qammod(bit_array(:, users_idx), M_QAM, 'InputType', 'bit'); 
% end
% 
% s_new = qammod(bit_array, M_QAM, 'InputType', 'bit');
% 
% disp(isequal(s_new, s_old));