% ------------------------------------------------------------------------------ %
% GENERAL PARAMETERS
% ------------------------------------------------------------------------------ %
M = 8;         % Number of base station antennas
K = 4;         % Number of users
N_BLK = 1;     % Number of blocks
B = 8;         % Number of bits
M_QAM = 2^B;   % Modulation (256 QAM) 
H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);  % Matriz de canal H MxK (8x4)

% ------------------------------------------------------------------------------ %
% MEMORY ALLOCATION
% ------------------------------------------------------------------------------ %
s = zeros(N_BLK, K);
bit_received = zeros(B*N_BLK, K);

% ------------------------------------------------------------------------------ %
% MODULATION AND PRE-PROCESSING
% ------------------------------------------------------------------------------ %

% Creating the random bit array
bit_array = randi([0,1], B*N_BLK, K);

% Creating the symbol matrix for each user N_BLKxK (1x4) 
for users_idx = 1:K
    s(:, users_idx) = qammod(bit_array(:,users_idx), M_QAM, 'InputType', 'bit');
end

% ------------------------------------------------------------------------------ %
% PRECODERS: ZF, MMSE, and MF
% ------------------------------------------------------------------------------ %

% ZF (Zero-Forcing) precoder
P_ZF = conj(H) / (H.' * conj(H));
x_ZF = P_ZF * s.';

% MMSE (Minimum Mean Square Error) precoder
sigma2 = 1; % Assume noise power (you can change this value)
P_MMSE = inv(H.' * H + sigma2 * eye(K)) * H.';
x_MMSE = P_MMSE * s.';

% MF (Matched Filter) precoder
P_MF = H.';
x_MF = P_MF * s.';

% Signal received
y_ZF = H.' * x_ZF;
y_MMSE = H.' * x_MMSE;
y_MF = H.' * x_MF;

% ------------------------------------------------------------------------------ %
% DEMODULATION AND BIT COMPARISON
% ------------------------------------------------------------------------------ %

% Demodulation for ZF
for users_idx = 1:K
    s_hat_ZF = y_ZF(users_idx, :).';
    bit_received(:, users_idx) = qamdemod(s_hat_ZF, M_QAM, 'OutputType', 'bit');
end

% Demodulation for MMSE
bit_received_MMSE = zeros(B*N_BLK, K); % Reset for MMSE
for users_idx = 1:K
    s_hat_MMSE = y_MMSE(users_idx, :).';
    bit_received_MMSE(:, users_idx) = qamdemod(s_hat_MMSE, M_QAM, 'OutputType', 'bit');
end

% Demodulation for MF
bit_received_MF = zeros(B*N_BLK, K); % Reset for MF
for users_idx = 1:K
    s_hat_MF = y_MF(users_idx, :).';
    bit_received_MF(:, users_idx) = qamdemod(s_hat_MF, M_QAM, 'OutputType', 'bit');
end

% ------------------------------------------------------------------------------ %
% BIT ERROR RATE (BER) CALCULATION
% ------------------------------------------------------------------------------ %

% Calculate BER for ZF
bit_errors_ZF = sum(bit_array ~= bit_received, 'all');
BER_ZF = bit_errors_ZF / numel(bit_array);

% Calculate BER for MMSE
bit_errors_MMSE = sum(bit_array ~= bit_received_MMSE, 'all');
BER_MMSE = bit_errors_MMSE / numel(bit_array);

% Calculate BER for MF
bit_errors_MF = sum(bit_array ~= bit_received_MF, 'all');
BER_MF = bit_errors_MF / numel(bit_array);

% Display the results
fprintf('ZF - Number of bit errors: %d, BER: %.5f\n', bit_errors_ZF, BER_ZF);
fprintf('MMSE - Number of bit errors: %d, BER: %.5f\n', bit_errors_MMSE, BER_MMSE);
fprintf('MF - Number of bit errors: %d, BER: %.5f\n', bit_errors_MF, BER_MF);
