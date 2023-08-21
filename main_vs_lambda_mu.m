clear all
clc;
close all
warning off

detector = [0,...% TS_ZF
    0,...% TS_MMSE
    0,...% TS_OSIC
    0,...% DL_TS_IS
    0,...% DL_TS_IS_AC
    0,...% DL_TS_IS_EM
    1];% DL_TS_IS_AC_EM

Nt = 32; Nr = Nt;
M_TS = 800; P_TS = round(M_TS/2);
n_iter = 1000;
line_style = {':r+', '-mo', '-k*', '-b*', '-ko', '-rs', '-ko'};
cutoff = 0.4; % reduce factor

% mu_vec = 1:1:5;
N = 2*Nt;
lambda_vec = [0:0.5:3].*1/N;
m_scheme = 'QPSK'; % BPSK, QPSK, 16QAM


global sim_ber sim_com
sim_ber = zeros(length(lambda_vec), length(detector));
sim_com = zeros(length(lambda_vec), length(detector));
for ss = 1:length(lambda_vec)
    snr_dB = 12;
    mu = 5;
    lambda = lambda_vec(ss)

    %disp(strcat('SNR = ', num2str(snr_dB), ' [dB]'));
    detect_mu(Nr, Nt, snr_dB, m_scheme, detector, n_iter, ss, M_TS, P_TS, cutoff, lambda, mu);
end
sim_ber;
sim_com;
plot_fig(lambda_vec, detector, m_scheme, Nt, Nr, 1, 0, 0, line_style)
% axis tight
hold all
