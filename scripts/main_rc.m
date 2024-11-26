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

BER = zeros(K, N_SNR);

% ####################################################################### %
% PARÂMETROS DE MODELAGEM DAS NÃO LINEARIDADES
% ####################################################################### %

A = 5;

% Amplificador de clipping ideal
clip = @(x, A0) min(abs(x), A0) .* exp(1j * angle(x));

% Amplificador valvulado de onda progressiva
kappa_A = 0.25;
chi_phi = 0.26;
kappa_phi = 0.25;

% Amplificador de estado sólido


% Polinomial


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
% TRANSMISSÃO
% ####################################################################### %

precoder_type = 'ZF'; % MF, MMSE

if strcmp(precoder_type, 'ZF')
    x = precoder_ZF * s.';

elseif strcmp(precoder_type, 'MF')
    x = precoder_MF * s.';

elseif strcmp(precoder_type, 'MMSE')
    x = zeros(M, N_BLK, N_SNR);
    for snr_idx = 1:N_SNR
        x(:,:,snr_idx) = precoder_MMSE(:,:,snr_idx) * s.';
    end

end

% Normalização de Potência de Transmissão
x_normalized = zeros(size(x));

if strcmp(precoder_type, 'MMSE')
    for m = 1:M
        for snr_idx = 1:N_SNR
            Px = norm(x(m, :, snr_idx))^2/size(x, 2);
            x_normalized(m, :, snr_idx) = x(m, :, snr_idx) / sqrt(Px);
        end
    end
else
    for m = 1:M
        Px = norm(x(m, :))^2/size(x, 2);
        x_normalized(m, :) = x(m, :) / sqrt(Px);
    end
end

% ####################################################################### %
% RECEPÇÃO
% ####################################################################### %

% non_linear_model = ideal, clipping, twta, solid amp, polinomial

v = sqrt(0.5) * (randn(K, N_BLK) + 1i*randn(K, N_BLK));
Pv = vecnorm(v,2,2).^2/N_BLK;
v_normalized = v./sqrt(Pv);

bit_received = zeros(B*N_BLK, K);

if strcmp(precoder_type, 'ideal')
    x = precoder_ZF * s.';

elseif strcmp(precoder_type, 'MF')
    x = precoder_MF * s.';

elseif strcmp(precoder_type, 'MMSE')
    x = zeros(M, N_BLK, N_SNR);
    for snr_idx = 1:N_SNR
        x(:,:,snr_idx) = precoder_MMSE(:,:,snr_idx) * s.';
    end

end

for snr_idx = 1:N_SNR

    if strcmp(precoder_type, 'ZF')
        y = H.' * sqrt(snr(snr_idx)) * x + v_normalized; % ideal
        % clipping
        % twta
        % solid amp
        % polinomial
        
    elseif strcmp(precoder_type, 'MF')
        y = H.' * sqrt(snr(snr_idx)) * x + v_normalized; % ideal
        % clipping
        % twta
        % solid amp
        % polinomial

    elseif strcmp(precoder_type, 'MMSE')
        y = H.' * (sqrt(snr(snr_idx)) * x(:,:,snr_idx)) + v_normalized; % ideal
        % clipping
        % twta
        % solid amp
        % polinomial
    end

    % Demodulação e cálculo da BER
    for users_idx = 1:K
        s_received = y(users_idx, :).';
        Ps_received = norm(s_received)^2/N_BLK;
        bit_received(:, users_idx) = qamdemod(sqrt(Ps(users_idx)/Ps_received) * s_received, M_QAM, 'OutputType', 'bit');
    end

    [~, BER(snr_idx)] = biterr(bit_received, bit_array);
end

% ####################################################################### %
% PLOTAGEM
% ####################################################################### %

% MÉDIA
figure;

BER_mean = mean(BER, 1);

set(gcf,'position',[0 0 800 600]);

switch precoder_type
    case 'ZF'
        semilogy(SNR, BER_mean,  'o-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(2,:)); 
        legend({'ZF'}, 'Location', 'northeast', 'FontSize', fontsize);
    case 'MF'
        semilogy(SNR, BER_mean,  's-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(3,:));
        legend({'MF'}, 'Location', 'northeast', 'FontSize', fontsize);
    case 'MMSE'
        semilogy(SNR, BER_mean,  'v-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(4,:));
        legend({'MMSE'}, 'Location', 'northeast', 'FontSize', fontsize);
end

xlabel('SNR (dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);

legend box off
set(gca, 'FontName', fontname, 'FontSize', fontsize);

% Salvando imagem
if savefig == 1
    saveas(gcf,'ber','fig');
    saveas(gcf,'ber','png');
    saveas(gcf,'ber','epsc2');
end