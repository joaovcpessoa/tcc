N_BLK = 1000;

M = 100;
K = 4;

B = 8;
M_QAM = 2^B;

SNR = -10:1:20;
N_SNR = length(SNR);

H = (randn(M, K) + 1i * randn(M, K)) / sqrt(2);

bit_array = randi([0,1], B*N_BLK, K);

% Resultado usando arrayfun
qam_modulation = @(bits) qammod(bits, M_QAM, 'InputType', 'bit'); 
s_cell = arrayfun(@(col_idx) qam_modulation(bit_array(:, col_idx)), 1:K, 'UniformOutput', false);
s_arrayfun = cell2mat(s_cell');
s_arrayfun = reshape(s_arrayfun, N_BLK, K);

% Resultado usando for loop
s_for = zeros(N_BLK, K);
for users_idx = 1:K
    s_for(:, users_idx) = qammod(bit_array(:,users_idx), M_QAM, 'InputType', 'bit');
end

% Comparar os resultados
disp(isequal(s_for, s_arrayfun))  % Deve retornar true se os resultados forem equivalentes