function [H, y, s, s_bin, s_NN_raw, sigma2, m_order] = ...
    gen_signal(Nt, Nr, m_scheme, snr, batch_S, batch_H, batch_Y, batch_S_NN, nn)

[symset, m_order] = loadAlphabet(m_scheme);
n_bit = log2(m_order);
Es = mean(abs(symset).^2); % % average symbol energy
sigma2 = Nt*Es / snr; % noise variance

s = batch_S(nn,:).';
y = batch_Y(nn,:).';
Htmp = batch_H(nn,:,:);
H = reshape(Htmp,2*Nr,2*Nt);
s_mod = s(1:Nt) + 1i*s(Nt+1:end);
s_int = qamdemod(s_mod,m_order);
s_bin = de2bi(s_int, n_bit);
s_NN_raw = squeeze(batch_S_NN(nn,:,:)).';
end