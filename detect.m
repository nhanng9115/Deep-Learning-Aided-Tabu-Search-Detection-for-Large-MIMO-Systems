function detect(Nr, Nt, snr_dB, m_scheme, detector, n_iter, ss, M, P, cutoff)

global sim_ber sim_com

N = 2*Nt;
if Nt == 16
    L = 10;
elseif Nt == 32
    L = 12;
elseif Nt == 32
    L = 15;
else
    L = 20;
end
iter_0 = 0;

doTS_ZF = detector(1);
doTS_MMSE = detector(2);
doTS_OSIC = detector(3);
doDLTS = detector(4); % initial solution + adaptive cutoff + efficient moves

ber_TS_ZF = zeros(n_iter-iter_0,1);
ber_TS_MMSE = zeros(n_iter-iter_0,1);
ber_TS_OSIC = zeros(n_iter-iter_0,1);
ber_DLTS = zeros(n_iter-iter_0,1);

com_TS_ZF = zeros(n_iter-iter_0,1);
com_TS_MMSE = zeros(n_iter-iter_0,1);
com_TS_OSIC = zeros(n_iter-iter_0,1);
com_DLTS = zeros(n_iter-iter_0,1);

%% Load data
snr = 10^(snr_dB/10);

if strcmp(m_scheme, 'QPSK')
    n_bit = 2;
    sig_file_name = strcat('.\test_data\QPSK\', num2str(Nt), '\', num2str(snr_dB),'dB.mat');
else
    n_bit = 4;
    sig_file_name = strcat('.\test_data\QAM\', num2str(Nt), '\', num2str(snr_dB),'dB.mat');
end

dat_signal = load(sig_file_name);

batch_S_NN = dat_signal.s_nn;
batch_S = dat_signal.s;
batch_Hb = dat_signal.H;
batch_Y = dat_signal.y;
parfor nn = 1:n_iter
    T = 0;
    ii = iter_0 + nn;
    [H, y, s, s_bin, s_NN_raw, sigma2, m_order] = ...
        gen_signal(Nt, Nr, m_scheme, snr, batch_S, batch_Hb, batch_Y, batch_S_NN, ii);
    
    if doTS_ZF% 'Conv. TS', use T and ET
        [s_lilear, n_oper_init] = init_detect(y, H, sigma2, m_scheme, 'ZF');
        [ber_TS_ZF(nn), com_TS_ZF(nn)] = ...
            TS_fast(y, H, M, P, m_scheme, s_lilear, m_order, n_bit, s_bin, nn, s, s_NN_raw, T, 1, cutoff, snr);
        com_TS_ZF(nn) = n_oper_init + com_TS_ZF(nn);
    end
    
    if doTS_MMSE% 'Conv. TS', use T and ET
        [s_lilear, n_oper_init] = init_detect(y, H, sigma2, m_scheme, 'MMSE');
        [ber_TS_MMSE(nn), com_TS_MMSE(nn)] = ...
            TS_fast(y, H, M, P, m_scheme, s_lilear, m_order, n_bit, s_bin, nn, s, s_NN_raw, T, 1, cutoff, snr);
        com_TS_MMSE(nn) = n_oper_init + com_TS_MMSE(nn);
    end
    
    if doTS_OSIC% 'Conv. TS', use T and ET
        [s_lilear, n_oper_init] = init_detect(y, H, sigma2, m_scheme, 'OSIC');
        [ber_TS_OSIC(nn), com_TS_OSIC(nn)] = ...
            TS_fast(y, H, M, P, m_scheme, s_lilear, m_order, n_bit, s_bin, nn, s, s_NN_raw, T, 1, cutoff, snr);
        com_TS_OSIC(nn) = n_oper_init + com_TS_OSIC(nn);
    end
    
    if doDLTS
        [ber_DLTS(nn), com_DLTS(nn)] = ...
            DL_TS(y, H, M, P, m_scheme, m_order, n_bit, s_bin, s, s_NN_raw, cutoff, snr, sigma2, 1, 1);
        n_oper_init = (2*N - 1 + 2*L)*N^2 + (2*N - 1 + 4*L)*N;
        com_DLTS(nn) = n_oper_init + com_DLTS(nn);
    end
end

% tmp
sim_ber(ss, :) = [mean(ber_TS_ZF), mean(ber_TS_MMSE), mean(ber_TS_OSIC), mean(ber_DLTS)]/(n_bit*Nt);
sim_com(ss, :) = [mean(com_TS_ZF), mean(com_TS_MMSE), mean(com_TS_OSIC), mean(com_DLTS)];

end % eof