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

N_BLK = 10000;

M = 100;
K = 4;
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);

B = 4;
M_QAM = 2^B;

SNR = -10:20;
N_SNR = length(SNR);
snr = 10.^(SNR/10);

% ####################################################################### %
% PARÂMETROS DE MODELAGEM DAS NÃO LINEARIDADES
% ####################################################################### %

% DESCREVER MELHOR ESSA VARIÁVEL
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
% MODELO DE TRANSMISSÃO (DOWNLINK)
% ####################################################################### %

s = zeros(N_BLK, K);
BER = zeros(K, N_SNR);
bit_array = randi([0,1], B*N_BLK, K);

for users_idx = 1:K
    s(:, users_idx) = qammod(bit_array(:,users_idx), M_QAM, 'InputType', 'bit'); 
end
     
Ps = vecnorm(s).^2/N_BLK;

precoder_ZF = conj(H) / (H.' * conj(H));
precoder_MF = conj(H) ./ (vecnorm(H).^2);

precoder_MMSE = zeros(M, K, N_SNR);
for snr_idx = 1:N_SNR
    precoder_MMSE(:,:,snr_idx)  = conj(H) / (H.' * conj(H) + 1/snr(snr_idx)*eye(K));
end

% IU
disp('Escolha o precodificador:');
disp('1 - ZF (Zero Forcing)');
disp('2 - MF (Matched Filter)');
disp('3 - MMSE (Minimum Mean Square Error)');

choice = input('Digite o número correspondente à sua escolha: ');

if choice == 1
    chosen_precoder = 'ZF';
    x = precoder_ZF * s.';
elseif choice == 2
    chosen_precoder = 'MF';
    x = precoder_MF * s.';
elseif choice == 3
    chosen_precoder = 'MMSE';
    x = zeros(M, N_BLK, N_SNR);
    for snr_idx = 1:N_SNR
        x(:,:,snr_idx) = precoder_MMSE(:,:,snr_idx) * s.';
    end
else
    error('Escolha inválida! Selecione 1, 2 ou 3.');
end

% Normalização de Potência de Transmissão
x_normalized = zeros(size(x));

if strcmp(chosen_precoder, 'MMSE')
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
% MODELO DE RECEPÇÃO (DOWNLINK)
% ####################################################################### %

v = sqrt(0.5) * (randn(K, N_BLK) + 1i*randn(K, N_BLK));
Pv = vecnorm(v,2,2).^2/N_BLK;
v_normalized = v./sqrt(Pv);

bit_received = zeros(B*N_BLK, K);

for snr_idx = 1:N_SNR
    if strcmp(chosen_precoder, 'ZF')
        y = H.' * sqrt(snr(snr_idx)) * x + v_normalized; 
    elseif strcmp(chosen_precoder, 'MF')
        y = H.' * sqrt(snr(snr_idx)) * x + v_normalized;
    elseif strcmp(chosen_precoder, 'MMSE')
        y = H.' * (sqrt(snr(snr_idx)) * x(:,:,snr_idx)) + v_normalized; 
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
% BER X SNR
% ####################################################################### %

% POR USUÁRIO
figure;

set(gcf,'position',[0 0 800 600]);

switch chosen_precoder
    case 'ZF'
        semilogy(SNR, BER,  'o-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(2,:)); 
        legend({'ZF'}, 'Location', 'northeast', 'FontSize', fontsize);
    case 'MF'
        semilogy(SNR, BER,  's-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(3,:));
        legend({'MF'}, 'Location', 'northeast', 'FontSize', fontsize);
    case 'MMSE'
        semilogy(SNR, BER,  'v-', 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(4,:));
        legend({'MMSE'}, 'Location', 'northeast', 'FontSize', fontsize);
end

xlabel('SNR (dB)', 'FontName', fontname, 'FontSize', fontsize);
ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);

legend box off
set(gca, 'FontName', fontname, 'FontSize', fontsize);

% Salvando imagem
if savefig == 1
    saveas(gcf,'ber_p_user','fig');
    saveas(gcf,'ber_p_user','png');
    saveas(gcf,'ber_p_user','epsc2');
end

% MÉDIA
figure;

BER_mean = mean(BER, 1);

set(gcf,'position',[0 0 800 600]);

switch chosen_precoder
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

% ####################################################################### %
% CONSTELAÇÃO DO SINAL
% ####################################################################### %

% Sinal decodificado (com normalização)
figure;

for users_idx = 1:K    
    s_received = y(users_idx, :).';
    Ps_received = norm(s_received)^2 / N_BLK;

    s_received_normalized = sqrt(Ps(users_idx) / Ps_received) * s_received;
    
    plot(real(s_received_normalized), imag(s_received_normalized),'.','MarkerSize', markersize,'Color',colors(users_idx, :));
    xlabel('Re', 'FontName', fontname, 'FontSize', fontsize);
    ylabel('Imag', 'FontName', fontname, 'FontSize', fontsize);
end