tic;
interm_folder_backup = interm_folder;
load(strcat(interm_folder,'/RQWfile.mat'));
RR = Rmatrices(:,:,targetCell, targetUser);

% ch_samples = 2;

dim = 2;

SII = eye(number_of_antennas)*(1/mu_val);
ZZ_2= zeros(number_of_antennas, ch_samples);
ZZ_3= zeros(ch_samples,1);
ZZ_4= zeros(number_of_antennas, number_of_antennas, ch_samples);
for nr = 1:length(NR_vec)
    clear W_Est;
    clear W_Est_diag;
    clear W_Est_diag_reg;
    
    rng shuffle;
    
    h_mat = generateh(R_sqrt_root,number_of_antennas,number_of_cells,number_of_users,ch_samples); %  (M, L, K, ch_samples)
    
    Ev = ZZ_2;
    Ev_diag = ZZ_2;
    Ev_diag_reg = ZZ_2;
    
    Eb = ZZ_2;
    Eb_diag = ZZ_2;
    Eb_diag_reg = ZZ_2;
    
    Evv = ZZ_3;
    Evv_diag = ZZ_3;
    Evv_diag_reg = ZZ_3;
    
    Evmat = ZZ_4;
    Evmat_diag = ZZ_4;
    Evmat_diag_reg = ZZ_4;
    
    Ebmat = ZZ_4;
    Ebmat_diag = ZZ_4;
    Ebmat_diag_reg = ZZ_4;
    trace_sum_mat = 0;
    
    for ii = 1:ch_samples
        sum_mat=0;
        sum_mat_diag=0;
        sum_mat_diag_reg=0;
        for k= 1:number_of_users
            load(strcat(interm_folder,'/RQest_',string(nr), '_', string(nq_sim),'_', string(k), '.mat'), 'W_Est', 'W_Est_diag','W_Est_diag_reg');
            [s1, h1] = mmse_cal(W_Est, TT, number_of_antennas, dim, pilotSequenceLength, mu_val, h_mat(:,:,k,ii));
            [s2, h2] = mmse_cal(W_Est_diag, TT, number_of_antennas, dim, pilotSequenceLength, mu_val, h_mat(:,:,k,ii));
            [s3, h3] = mmse_cal(W_Est_diag_reg, TT, number_of_antennas, dim, pilotSequenceLength, mu_val, h_mat(:,:,k,ii));
            clear W_Est;
            clear W_Est_diag;
            clear W_Est_diag_reg;
            if k== targetUser
                h_target = h1;
                h_target_diag = h2;
                h_target_diag_reg = h3;
            end
            sum_mat = sum_mat + s1;
            sum_mat_diag = sum_mat_diag + s2;
            sum_mat_diag_reg = sum_mat_diag_reg + s3;
        end
        trace_sum_mat = trace_sum_mat + (1/ch_samples)*sum(sum_mat(:));
        
        for tt = 1:TT 
            v = compute_v_or_b(sum_mat(:,:,tt), h_target(:,tt), SII, 0); % M, 1
            v_diag = compute_v_or_b(sum_mat_diag(:,:,tt), h_target_diag(:,tt), SII, 0); % M, 1
            v_diag_reg = compute_v_or_b(sum_mat_diag_reg(:,:,tt), h_target_diag_reg(:,tt), SII, 0); % M, 1
            
            Ev(:,ii) = Ev(:,ii) + (1/TT)*v;
            Ev_diag(:,ii) = Ev_diag(:,ii) + (1/TT)*v_diag;
            Ev_diag_reg(:,ii) = Ev_diag_reg(:,ii) + (1/TT)*v_diag_reg;
            
            Evv(ii) = Evv(ii) + (1/TT)*(v'*v);
            Evv_diag(ii) = Evv_diag(ii) + (1/TT)*(v_diag'*v_diag);
            Evv_diag_reg(ii) = Evv_diag_reg(ii) + (1/TT)*(v_diag_reg'*v_diag_reg);
            
            b = compute_v_or_b(sum_mat(:,:,tt), h_target(:,tt), SII, 1);
            b_diag = compute_v_or_b(sum_mat_diag(:,:,tt), h_target_diag(:,tt), SII, 1);
            b_diag_reg = compute_v_or_b(sum_mat_diag_reg(:,:,tt), h_target_diag_reg(:,tt), SII, 1);
            
            Eb(:,ii) = Eb(:,ii) + (1/TT)*b;
            Eb_diag(:,ii) = Eb_diag(:,ii) + (1/TT)*b_diag;
            Eb_diag_reg(:,ii) = Eb_diag_reg(:,ii) + (1/TT)*b_diag_reg;
            
            Evmat(:, :, ii) = Evmat(:, :, ii) + (1/TT)*(v*v');
            Evmat_diag(:, :, ii) = Evmat_diag(:, :, ii) + (1/TT)*(v_diag*v_diag');
            Evmat_diag_reg(:, :, ii) = Evmat_diag_reg(:, :, ii) + (1/TT)*(v_diag_reg*v_diag_reg');
            
            Ebmat(:, :, ii) = Ebmat(:, :, ii) + (1/TT)*(b*b');
            Ebmat_diag(:, :, ii) = Ebmat_diag(:, :, ii) + (1/TT)*(b_diag*b_diag');
            Ebmat_diag_reg(:, :, ii) = Ebmat_diag_reg(:, :, ii) + (1/TT)*(b_diag_reg*b_diag_reg');
        end
    end
    
    
    NR = NR_vec(nr);

    num_hat = 0;
    num_hat_diag = 0;
    num_hat_diag_reg = 0;
    
    num_hat_DL = 0;
    num_hat_diag_DL = 0;
    num_hat_diag_reg_DL = 0;

    den1_hat = 0;
    den1_hat_diag = 0;
    den1_hat_diag_reg = 0;

    den1_hat_DL = 0;
    den1_hat_diag_DL = 0;
    den1_hat_diag_reg_DL = 0;

    den2_hat = 0;
    den2_hat_diag = 0;
    den2_hat_diag_reg = 0;
    
    for ii = 1:ch_samples
        num_hat = num_hat + (1/ch_samples) * (Ev(:,ii)'*h_mat(:,targetCell,targetUser, ii));
        num_hat_diag = num_hat_diag + (1/ch_samples) * (Ev_diag(:,ii)'*h_mat(:,targetCell,targetUser, ii));
        num_hat_diag_reg = num_hat_diag_reg + (1/ch_samples) * (Ev_diag_reg(:,ii)'*h_mat(:,targetCell,targetUser, ii));

        num_hat_DL = num_hat_DL + (1/ch_samples) * (Eb(:,ii)'*h_mat(:,targetCell,targetUser, ii));
        num_hat_diag_DL = num_hat_diag_DL + (1/ch_samples) * (Eb_diag(:,ii)'*h_mat(:,targetCell,targetUser, ii));
        num_hat_diag_reg_DL = num_hat_diag_reg_DL + (1/ch_samples) * (Eb_diag_reg(:,ii)'*h_mat(:,targetCell,targetUser, ii));
        
        hh = zeros(number_of_antennas, number_of_antennas);

        for l=1:number_of_cells
            for k=1:number_of_users
                hh = hh + (h_mat(:,l,k, ii)*h_mat(:,l,k, ii)');
            end
        end
        
        den1_hat = den1_hat + (1/ch_samples) * trace(hh * Evmat(:, :, ii));
        den1_hat_diag = den1_hat_diag + (1/ch_samples) * trace(hh * Evmat_diag(:, :, ii));
        den1_hat_diag_reg = den1_hat_diag_reg + (1/ch_samples) * trace(hh * Evmat_diag_reg(:, :, ii));

        den1_hat_DL = den1_hat_DL + (1/ch_samples) * trace(hh * Ebmat(:, :, ii));
        den1_hat_diag_DL = den1_hat_diag_DL + (1/ch_samples) * trace(hh * Ebmat_diag(:, :, ii));
        den1_hat_diag_reg_DL = den1_hat_diag_reg_DL + (1/ch_samples) * trace(hh * Ebmat_diag_reg(:, :, ii));
        
        den2_hat = den2_hat + (1/ch_samples) * (1/mu_val) * (Evv(ii));
        den2_hat_diag = den2_hat_diag + (1/ch_samples) * (1/mu_val) * (Evv_diag(ii));
        den2_hat_diag_reg = den2_hat_diag_reg + (1/ch_samples) * (1/mu_val) * (Evv_diag_reg(ii));
    end
    if nr==5 && nq_sim==6
        trace_Ev = sum(mean(Ev, 2));
        if abs(den1_hat)>10
            disp('--------------------');
            disp('block_ID', string(block_ID));
            disp('trace_sum_mat', string(trace_sum_mat));
            disp('trace_Ev', string(trace_Ev));
            disp('den1_hat', string(sum(den1_hat(:))));
        elseif abs(den1_hat)>0.9 && abs(den1_hat)< 2
            disp('--------------------');
            disp('block_ID', string(block_ID));
            disp('trace_sum_mat', string(trace_sum_mat));
            disp('trace_Ev', string(trace_Ev));
            disp('den1_hat', string(sum(den1_hat(:))));
        end
    end

    
    save(strcat(out_folder,'/outs_ZF_',string(nr), '_', string(nq_sim), '_',  string(block_ID), '.mat'), 'num_hat', 'num_hat_diag',...
        'num_hat_diag_reg', 'num_hat_DL', 'num_hat_diag_DL', 'num_hat_diag_reg_DL', 'den1_hat', 'den1_hat_diag', 'den1_hat_diag_reg', ...
        'den1_hat_DL', 'den1_hat_diag_DL', 'den1_hat_diag_reg_DL', 'den2_hat', 'den2_hat_diag', 'den2_hat_diag_reg');


    ['ZF' string(nr) string(toc)]
end

function v_or_b = compute_v_or_b(pre_mat, h_target, SII, ULDL_flag)
    v = (pre_mat+SII)\h_target;
%     v = h_target;
    const = 1; % TODO
    if ULDL_flag == 0
        v_or_b = v;
    else
        v_or_b = v/const;
    end
end

function [sum_mat, h_out]= mmse_cal(W_Est, TT, number_of_antennas, dim, P, mu, h_mat)
%     h_mat is of size (M,L)
    rng shuffle;
    ns = crandn(number_of_antennas,TT);

    h_sum = sum(h_mat, dim);
    
    for tt = 1:TT
        W_est = W_Est(:,:,tt);
        h_LS = h_sum + (1/sqrt(P*mu)) * ns(:,tt);

        h_mmse = W_est*h_LS;
        sum_mat(:,:,tt) = h_mmse*h_mmse'; % sum over k
        h_out(:,tt) = h_mmse;

    end
end