clear all
clc;
% close all
warning off

detector = [1,...% TS_ZF
    1,...% TS_MMSE
    1,...% TS_OSIC
    1];% DL_TS

test_mode = 1;
if test_mode == 1 % QPSK
    SNR_dB_vec = 12;
    m_scheme = 'QPSK'; % BPSK, QPSK, 16QAM
    Nt = 32; Nr = Nt;
    M_TS = 800; P_TS = round(M_TS/2);
    n_iter_vec = 100;
else % 16-QAM ===================================================
    SNR_dB_vec = 10:2:22;
    m_scheme = '16QAM'; % BPSK, QPSK, 16QAM
    Nt = 8; Nr = Nt;
    M_TS = 2500; P_TS = round(M_TS/2);
    n_iter_vec = 4000.+[200,400,800,1000,2000,4000,6000];
end

% line_style = {'-mo', ':r+', '-k*', '-b*', '-ko', '-rs', '--ro'};
line_style = {':r+', '-mo', '-k*', '-k*', '-ko', '-rs', '-ko'};

cutoff = 0.4; % reduce factor

global sim_ber sim_com sim_dis 
sim_ber = zeros(length(SNR_dB_vec), length(detector));
sim_com = zeros(length(SNR_dB_vec), length(detector));
sim_dis = zeros(length(SNR_dB_vec), length(detector));
for ss = 1:length(SNR_dB_vec)
    snr_dB = SNR_dB_vec(ss);
    disp(strcat('SNR = ', num2str(snr_dB), ' [dB]'));
    detect_time(Nr, Nt, snr_dB, m_scheme, detector, n_iter_vec(ss), ss, M_TS, P_TS, cutoff);
end
sim_ber
sim_com
plot_fig(SNR_dB_vec, detector, m_scheme, Nt, Nr, 1, 1, 0, line_style)
% axis tight
hold all
