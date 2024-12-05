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
root_save = ['C:\Users\joaov_zm1q2wh\OneDrive\Code\github\tcc\scripts\images\'];

% ####################################################################### %
%% PARAMETROS PRINCIPAIS
% ####################################################################### %

load('ber_mf.mat');

precoder_type = 'MF';
amplifiers_type = {'IDEAL', 'CLIP', 'TWT', 'SS'};
A0 = [0.5, 1.0, 1.5, 2.0, 2.5];

% ####################################################################### %
%% REGRESSAO SIMBOLICA
% ####################################################################### %

snr_linear = 10 .^ (SNR / 10);
ber = mean(BER(:, :, 1, 1), 1);

syms a b c d x

f = a * exp(-b * x) + c * x ^ d;

x_data = snr_linear;
y_data = ber;
params = fitnlm(x_data, y_data, @(b, x) b(1) * exp(-b(2) * x) + b(3) * x .^ b(4), [1, 1, 1, 1]);

coeff = params.Coefficients.Estimate;

x_fit = linspace(min(x_data), max(x_data), 1000);
y_fit = coeff(1) * exp(-coeff(2) * x_fit) + coeff(3) * x_fit .^ coeff(4);

disp(x_data);
disp(y_data);
disp(params);
disp(x_fit);
disp(y_fit);

% ####################################################################### %
%% PLOT
% ####################################################################### %

for amp_idx = 1:N_AMP
    figure;
    set(gcf, 'position', [0 0 800 600]);
    
    for a_idx = 1:N_A0
        semilogy(SNR, mean(BER(:,:,amp_idx,a_idx),1), 'LineWidth', linewidth, 'MarkerSize', markersize, 'Color', colors(a_idx+1,:));
        hold on;
    end

    xlabel('SNR (dB)', 'FontName', fontname, 'FontSize', fontsize);
    ylabel('BER', 'FontName', fontname, 'FontSize', fontsize);
    title(sprintf('Amplificador: %s', amplifiers_type{amp_idx}), 'FontName', fontname, 'FontSize', fontsize);
    
    legend(arrayfun(@(a) sprintf('A=%.1f', a), A0, 'UniformOutput', false), ...
           'Location', 'northeast', 'FontSize', fontsize);
    legend box off;
    
    set(gca, 'FontName', fontname, 'FontSize', fontsize);

    graph_name = sprintf('BER_%s_%s.fig', precoder_type, amplifiers_type{amp_idx});
    
    if savefig == 1
        saveas(gcf,[root_save graph_name],'fig');
        saveas(gcf,[root_save graph_name],'png');
        saveas(gcf,[root_save graph_name],'epsc2');
    end
end
