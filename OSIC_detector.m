function [X_hat, n_add, n_mul]=OSIC_detector(y,H,sigma2,N,OSIC_type,modScheme)

% MIMO-OFDM Wireless Communications with MATLAB,  Yong Soo Cho
% Edited by Nguyen Thanh Nhan
X_hat = zeros(N,1);
n_add = 0;
n_mul = 0;
switch OSIC_type
    case 1  % Post_detection_SINR
        Order=[];  % detection order
        index_array = 1:N; % yet to be detected signal index
        % complexity to compute H'*H + sigma2*eye(k) in advance
        n_add = n_add + 3*N^3 - N^2 - 2*N + 2;
        n_mul = n_add + 3*N^3 + N^2 - 2*N;
        % V-BLAST
        for stage = 1:N
            k = N-(stage-1);
            Wmmse = (H'*H + sigma2*eye(size(H,2)))^(-1) * H';  % MMSE filter
            % complexity for matrix inversion
            n_mul = n_mul + 1/3*k^3 + 1/2*k^2 + 1/6*k;
            n_add = n_add + 1/3*k^3 + 1/2*k^2 + 1/6*k - 1;
        
            WmmseH = Wmmse*H;
            n_add = n_add + k^3 - k^2;
            n_mul = n_mul + k^3;
        
            SINR=[];
            
            for i = 1:k
                P_S = abs(WmmseH(i,i))^2;
                P_IN = norm(WmmseH(i,[1:i-1, i+1:k]))^2 + sigma2*norm(Wmmse(i,:))^2;
                SINR(i) = P_S/P_IN; % SINR calculation
                
                n_mul = n_mul + 2*k + 4;
                n_add = n_add + 2*k;
            end
            [SINR_max, index_temp] = max(SINR);    % ordering using SINR
            Order = [Order, index_array(index_temp)];
            index_array = index_array([1:index_temp-1 index_temp+1:end]);
            x_temp(stage) = Wmmse(index_temp,:)*y;     % MMSE filtering
            n_mul = n_mul + N;
            n_add = n_add + N;
            
            X_hat(stage) = mod_slicer(x_temp(stage), modScheme);
            y_tilde = y - H(:,index_temp)*X_hat(stage); % interference subtraction
            n_mul = n_mul + N;
            n_add = n_add + N;
                
            H_tilde = H(:,[1:index_temp-1 index_temp+1:k]); % new H
            H = H_tilde;   
            y = y_tilde;
        end
        X_hat(Order) = X_hat;
    case 2 % column_norm ordering detection
        X_hat = zeros(N,1);
        G = inv(H);           % inverse of H
        n_mul = n_mul + 1/3*N^3 + 1/2*N^2 + 1/6*N;
        n_add = n_add + 1/3*N^3 + 1/2*N^2 + 1/6*N - 1;
        for i=1:N            % column_norm calculation
            norm_array(i) = norm(H(:,i));
            n_mul = n_mul + N;
            n_add = n_add + N-1;
        end
        [sorted_norm_array,Order_temp] = sort(norm_array);
        Order = fliplr(Order_temp);
        % V-BLAST
        for stage=1:N
            x_temp=G(Order(stage),:)*y;    % Tx signal estimation
%             stage
%             size(x_temp)
            X_hat(stage) = mod_slicer(x_temp, modScheme);
            y_tilde = y - H(:,Order(stage))*X_hat(Order(stage));
            n_mul = n_mul + N;
            n_add = n_add + N;
        end
    otherwise % OSIC with Post_detection_SNR ordering
        Order = [];
        index_array = 1:N; % set of indices of signals to be detected
        % V-BLAST
        for stage=1:N
            k = N-(stage-1);
            
            W_ZF = (H'*H)^(-1)*H'; % ZF combining matrix
            n_mul = n_mul + 1/3*k^3 + 1/2*k^2 + 1/6*k;
            n_add = n_add + 1/3*k^3 + 1/2*k^2 + 1/6*k - 1;
            norm_array = [];
            
            for i=1:k % detection ordering
                norm_array(i) = norm(W_ZF(i,:));
                n_mul = n_mul + N;
                n_add = n_add + N-1;
            end
            [val_min,index_min] = min(norm_array); % ordering in SNR
            Order = [Order index_array(index_min)];
            index_array = index_array([1:index_min-1 index_min+1:end]);
            x_temp(stage) = W_ZF(index_min,:)*y;  % Tx signal estimation
            n_mul = n_mul + N;
            n_add = n_add + N-1;
            X_hat(stage) = mod_slicer(x_temp(stage), modScheme);
            y_tilde = y - H(:,index_min)*X_hat(stage); % interference subtraction
            n_mul = n_mul + N;
            n_add = n_add + N;
            H_tilde = H(:,[1:index_min-1 index_min+1:N-(stage-1)]); % new H
            H = H_tilde;   y = y_tilde;
        end
        X_hat(Order) = X_hat;
end
% X_hat = X_hat.';