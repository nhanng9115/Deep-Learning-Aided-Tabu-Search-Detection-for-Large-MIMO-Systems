function [ber, comp, dis] = ...
    TS_fast(y, H, M, P, m_scheme, s_init, m_order, n_bit, x_bin, nn, s, s_NN_raw, T, use_T, cutoff, snr)

N = size(H,2);
K = size(H,1);
Nt = N/2;
% counting number of complex additions, multiplications, comparisons
n_mul = 0; n_add = 0; n_comp = 0;

m_init = norm(y - H*s_init)^2; % = disCand'*disCand => need Nr muls and Nr-1 adds

cand = s_init;
cand_c = cand(1:Nt) + 1i*cand(Nt+1:end);
tb_list = cand_c; % Put the initial solution to tabu list
tb_metric = m_init;

iter_max = round(M*cutoff);
s_NN = mod_slicer(s_NN_raw, m_scheme).';%sign(s_NN_raw);
err_threshold = min(1/snr, 0.5);
err = abs(s_NN - s_NN_raw);
err_pos = find(err > err_threshold);
n_err = length(err_pos);

%% test
s_hat = cand;
m_hat = m_init;
p = 1;
count = 0;
iter = 1;
m_best = [];
dis = 0;
% s - s_init

while iter <= M && count <= iter_max
    
    %% Find neighbor set
    nb_set = find_neighbors(cand_c, m_scheme,tb_list, Nt);
    L = size(nb_set, 2);
    if L == 0
        break
    end

    %% Compute neighbors metrics
    list_m = vecnorm(y - H*nb_set).^2;
    % counting
    n_mul = n_mul + L*2*K;
    n_add = n_add + L*2*K-1;
    
    %% Find the best neighbor
    [m_cand, i_min] = min(list_m); % L comps
    n_comp = n_comp + L;
    cand = nb_set(:,i_min);


    %% save the best solution
    p = p + 1; % increase the tabu lenth
    if m_cand < m_hat
        s_hat = cand;
%         (s-s_hat).'
        m_hat = m_cand;
        count = 0;
    else
        count = count + 1;
    end
    m_best = cat(2,m_best,m_cand);
    
    if iter == n_err
        dis = dis + norm(s_hat - s)^2;
    end
    %% release the first solution if tabulist is full and push a new solution
    if p > P
        tb_list = tb_list(:,2:end);
        tb_metric = tb_metric(:,2:end);
        p = p - 1;
    end
    cand_c = cand(1:Nt) + 1i*cand(Nt+1:end);
    tb_list = cat(2, tb_list, cand_c);
    tb_metric = cat(2, tb_metric, m_cand);
    
    iter = iter + 1;
end
% (s - s_init).'
% (s - s_hat).'
% if ~isequal(s, s_hat)
%     disp(strcat('#sample = ', num2str(nn)));
%     len = length(find(s - s_hat));
%     disp(strcat('#error symbols = ', num2str(len)));
% end
% disp('======================================================================================');

s_hat_c = s_hat(1:Nt) + 1i*s_hat(Nt+1:end);
x_hat_demod = qamdemod(s_hat_c,m_order);
x_hat_bin = de2bi(x_hat_demod, n_bit);
ber = biterr(x_bin,x_hat_bin); % ber

comp = n_mul + n_add; % complexity
end % end function