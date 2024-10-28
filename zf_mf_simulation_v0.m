
% ------------------------------------------------------------------------------ %
% GENERAL PARAMETERS
% ------------------------------------------------------------------------------ %
M = 500;         % Number of base station antennas
K = 4;         % Number of users
N_BLK = 1;     % Number of blocks
B = 8;         % Number of bits
M_QAM = 2^B;   % Modulation (256 QAM) 
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);  % Matriz de canal H MxK (8x4)

% ------------------------------------------------------------------------------ %
% MEMORY ALLOCATION
% ------------------------------------------------------------------------------ %
s = zeros(N_BLK, K);

% ------------------------------------------------------------------------------ %
% MODULATION AND PRE-PROCESSING
% ------------------------------------------------------------------------------ %

% Creating the random bit array
bit_array = randi([0,1], B*N_BLK, K);

% Creating the symbol matrix for each user N_BLKxK (1x4) 
for users_idx = 1:K
    s(:, users_idx) = qammod(bit_array(:,users_idx), M_QAM, 'InputType', 'bit');
end
     
P_ZF = conj(H) / (H.' * conj(H));  % ZF (Zero-forcing) precoder
P_MF = conj(H);                    % MF (Matched Filter) precoder

% Signal transmitted
x_ZF = P_ZF * s.';
x_MF = P_MF * s.';

% Signal received
y_ZF = H.'* x_ZF;
y_MF = H.'* x_MF;

% ------------------------------------------------------------------------------ %
% ADDING AWGN NOISE
% ------------------------------------------------------------------------------ %

% SNR_dB = 20;             % SNR (dB)
% SNR = 10^(SNR_dB / 10);  % Linear scale
% 
% % Power of transmitted signal = norm of transmitted signal
% % The noise power is then determined based on the SNR
% 
% % ZF
% signal_power_ZF = norm(x_ZF, 'fro')^2 / (M*K);  % Average power of transmitted signal
% noise_power_ZF = signal_power_ZF / SNR;         % Noise power for ZF
% 
% % MF
% signal_power_MF = norm(x_MF, 'fro')^2 / (M*K);  % Average power of transmitted signal
% noise_power_MF = signal_power_MF / SNR;         % Noise power for MF
% 
% % Generate AWGN noise
% noise_ZF = sqrt(noise_power_ZF) * (randn(M, 1) + 1i * randn(M, 1)) / sqrt(2);
% noise_MF = sqrt(noise_power_MF) * (randn(M, 1) + 1i * randn(M, 1)) / sqrt(2);
% 
% % Add noise to received signals
% y_ZF_with_noise = y_ZF + noise_ZF;  % Received signal with noise for ZF
% y_MF_with_noise = y_MF + noise_MF;  % Received signal with noise for MF

% AWGN com variância 10^-3


% ------------------------------------------------------------------------------ %
% DEMODULATION AND BIT COMPARISON
% ------------------------------------------------------------------------------ %
% ZF
bit_received_ZF = zeros(B*N_BLK, K);
for users_idx = 1:K
    s_received_ZF = y_ZF(users_idx, :).';
    bit_received_ZF(:, users_idx) = qamdemod(s_received_ZF, M_QAM, 'OutputType', 'bit');
end

% MF
bit_received_MF = zeros(B*N_BLK, K); % Reset for MF
for users_idx = 1:K
    s_received_MF = y_MF(users_idx, :).';
    bit_received_MF(:, users_idx) = qamdemod(s_received_MF, M_QAM, 'OutputType', 'bit');
end

% ------------------------------------------------------------------------------ %
% BIT ERROR RATE (BER) CALCULATION
% ------------------------------------------------------------------------------ %

noise_variance = 10^-3;

noise_ZF = sqrt(noise_variance/2) * (randn(size(s_received_ZF)) + 1i * randn(size(s_received_ZF)));
noise_MF = sqrt(noise_variance/2) * (randn(size(s_received_MF)) + 1i * randn(size(s_received_MF)));

received_signal_with_noise_ZF = s_received_ZF + noise_ZF;
received_signal_with_noise_MF = s_received_MF + noise_MF;

% BER for ZF
bit_errors_ZF = sum(bit_array ~= bit_received, 'all');
BER_ZF = bit_errors_ZF / numel(bit_array);

% BER for MF
bit_errors_MF = sum(bit_array ~= bit_received_MF, 'all');
BER_MF = bit_errors_MF / numel(bit_array);

% Results
disp('ZF:')
disp('Sinal recebido através do canal H:')
disp(s_received_ZF);
disp('Sinal recebido através do canal H (com ruído):');
disp(received_signal_with_noise_ZF);
fprintf('ZF - Number of bit errors: %d, BER: %.5f\n', bit_errors_ZF, BER_ZF);

disp('')
disp('--------------------------------------------')
disp('')

disp('MF:')
disp('Sinal recebido através do canal H:')
disp(s_received_ZF);
disp('Sinal recebido através do canal H (com ruído):');
disp(received_signal_with_noise_MF);
fprintf('MF - Number of bit errors: %d, BER: %.5f\n', bit_errors_MF, BER_MF);