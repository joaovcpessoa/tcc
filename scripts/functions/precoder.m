function [y] = precoder(x, type)

precoder_ZF = conj(H) / (H.' * conj(H));
precoder_MF = conj(H) ./ (vecnorm(H).^2);
precoder_MMSE = zeros(M, K, N_SNR);
for snr_idx = 1:N_SNR
    precoder_MMSE(:,:,snr_idx) = conj(H) / (H.' * conj(H) + 1/snr(snr_idx) * eye(K));
end

        y_ZF = H.' * amplifiers{pre_idx}(sqrt(snr(snr_idx)) * precoder_ZF * s.') + v_normalized; 
        y_MF = H.' * amplifiers{pre_idx}(sqrt(snr(snr_idx)) * precoder_MF * s.') + v_normalized;
        y_MMSE = H.' * amplifiers{pre_idx}(sqrt(snr(snr_idx)) * precoder_MMSE(:,:,snr_idx) * s.') + v_normalized; 