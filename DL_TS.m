function [ber, comp] = ...
    DL_TS(y, H, M, P, m_scheme, m_order, n_bit, x_bin, s, s_NN_raw, cutoff, snr, sigma2, use_AC, use_EM)

N = size(H,2);
K = size(H,1);
Nt = N/2;
n_mul = 0; n_add = 0;

%% NN-based initial solution ==================================
s_NN = mod_slicer(s_NN_raw, m_scheme).';
s_init = s_NN;

err = abs(s_NN - s_NN_raw);
if strcmp(m_scheme,'QPSK') && Nt == 16
    mu = 4;
elseif strcmp(m_scheme,'QPSK') && Nt == 24
    mu = 10;
elseif strcmp(m_scheme,'QPSK') && Nt == 32
    mu = 5;
else
    mu = 5;
end
err_threshold = min(sigma2/N, 0.5);


err_pos = find(err > err_threshold);
n_err = length(err_pos);
iter_max = M * min(cutoff, mu*n_err/N) + 2;


%% TS initialization with init solution
cand = s_init;
cand_c = cand(1:Nt) + 1i*cand(Nt+1:end);
tb_list = cand_c;
tb_metric = norm(y - H*cand)^2;

%% To save and update the final solution
s_hat = cand;
m_hat = tb_metric;
count = 0;
iter = 1;
considered_pos = [];
while iter <= M && count < iter_max
    
    % Find nb set
    nb_set = find_neighbors(cand_c, m_scheme, tb_list, Nt) ;
    L = size(nb_set, 2);
    if L == 0
        break
    end
    
    if (use_EM == 1) && (iter <= 4*n_err) %&& (count <= 2)
        
        list_m = Inf*ones(1,L);
        for ii = 1:L
            x = nb_set(:,ii);
            d = find(x - cand);
            if ismember(d,err_pos)
                m_x = norm(y - H*x)^2;
                list_m(ii) = m_x;
                
                % counting
                n_mul = n_mul + 2*K;
                n_add = n_add + 2*K-1;
            end
        end
        [m_cand, i_min] = min(list_m);
        best_nb = nb_set(:,i_min);
        d_min = find(best_nb - cand);
        if ismember(d_min,considered_pos)
            %iter = iter - 1;
        end
        
        considered_pos = cat(2,considered_pos,d_min);
        cand = best_nb;
        
    else
        list_m = vecnorm(y - H*nb_set).^2;
        n_mul = n_mul + L*2*K;
        n_add = n_add + L*2*K-1;
        
        [m_cand, i_min] = min(list_m);
        best_nb = nb_set(:,i_min);
        cand = best_nb;
    end
    
    % update the best solution
    if m_cand < m_hat
        s_hat = cand;
        m_hat = m_cand;
        count  = 0;
    else
        count = count + 1;
    end
    %m_hat
    
    %% release the first solution if tabulist is full and push a new solution
    if length(tb_metric) > P
        tb_list = tb_list(:,2:end);
        tb_metric = tb_metric(:,2:end);
    end
    
    cand_c = cand(1:Nt) + 1i*cand(Nt+1:end);
    tb_list = cat(2, tb_list, cand_c);
    tb_metric = cat(2, tb_metric, m_cand);
    iter = iter + 1;
end

s_hat_c = s_hat(1:Nt) + 1i*s_hat(Nt+1:end);
x_hat_demod = qamdemod(s_hat_c,m_order);
x_hat_bin = de2bi(x_hat_demod, n_bit);
ber = biterr(x_bin,x_hat_bin);

comp = n_mul + n_add;

end % end function