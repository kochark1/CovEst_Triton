tic;
load(strcat(interm_folder,'/RQWfile.mat'));
RR = Rmatrices(:,:,targetCell, targetUser);
for nr = 1:length(NR_vec)
    clear W_Est;
    clear W_Est_diag;
    NR = NR_vec(nr);

    load(strcat(interm_folder,'/RQest_',string(nr), '_', string(nq_sim),'_', string(targetUser), '.mat'), 'W_Est', 'W_Est_diag','W_Est_diag_reg');


    num_hat_t = zeros(TT,1);
    num_hat_diag_t = zeros(TT,1);
    num_hat_diag_reg_t = zeros(TT,1);

    den1_hat_t = zeros(TT,1);
    den1_hat_diag_t = zeros(TT,1);
    den1_hat_diag_reg_t = zeros(TT,1);

    den1_hat_DL_t = zeros(TT,1);
    den1_hat_diag_DL_t = zeros(TT,1);
    den1_hat_diag_reg_DL_t = zeros(TT,1);

    den2_hat_t = zeros(TT,1);
    den2_hat_diag_t = zeros(TT,1);
    den2_hat_diag_reg_t = zeros(TT,1);

    for t = 1:TT
        W_est = W_Est(:,:,t);
        W_est_diag = W_Est_diag(:,:,t);
        W_est_diag_reg = W_Est_diag_reg(:,:,t);

        num_hat_t(t) = trace(W_est'*RR);
        num_hat_diag_t(t) = trace(W_est_diag'*RR);
        num_hat_diag_reg_t(t) = trace(W_est_diag_reg'*RR);

        den1_hat_t(t) = trace(W_est*Qmatrix*W_est'*Rsum);
        den1_hat_diag_t(t) = trace(W_est_diag*Qmatrix*W_est_diag'*Rsum);
        den1_hat_diag_reg_t(t) = trace(W_est_diag_reg*Qmatrix*W_est_diag_reg'*Rsum);

        den1_hat_DL_t(t) = trace(W_est*Qmatrix*W_est'*Rsum_DL);
        den1_hat_diag_DL_t(t) = trace(W_est_diag*Qmatrix*W_est_diag'*Rsum_DL);
        den1_hat_diag_reg_DL_t(t) = trace(W_est_diag_reg*Qmatrix*W_est_diag_reg'*Rsum_DL);
        for l=1:number_of_cells
            den2_hat_t(t) = den2_hat_t(t) + abs(trace(W_est'*Rmatrices(:,:,l, targetUser)))^2;
            den2_hat_diag_t(t) = den2_hat_diag_t(t) + abs(trace(W_est_diag'*Rmatrices(:,:,l, targetUser)))^2;
            den2_hat_diag_reg_t(t) = den2_hat_diag_reg_t(t) + abs(trace(W_est_diag_reg'*Rmatrices(:,:,l, targetUser)))^2;
        end
    end

    num_hat = abs(mean(num_hat_t))^2;
    num_hat_diag = abs(mean(num_hat_diag_t))^2;
    num_hat_diag_reg = abs(mean(num_hat_diag_reg_t))^2;

    den_hat = mean(den1_hat_t + den2_hat_t) - num_hat;
    den_hat_diag = mean(den1_hat_diag_t + den2_hat_diag_t) - num_hat_diag;
    den_hat_diag_reg = mean(den1_hat_diag_reg_t + den2_hat_diag_reg_t) - num_hat_diag_reg;

    den_hat_DL = mean(den1_hat_DL_t + den2_hat_t) - num_hat  +(1/lambda_DL);
    den_hat_diag_DL = mean(den1_hat_diag_DL_t + den2_hat_diag_t) - num_hat_diag +(1/lambda_DL);
    den_hat_diag_reg_DL = mean(den1_hat_diag_reg_DL_t + den2_hat_diag_reg_t) - num_hat_diag_reg +(1/lambda_DL);

    save(strcat(out_folder,'/outs_',string(nr), '_', string(nq_sim), string(block_ID),'.mat'), 'num_hat', 'num_hat_diag','num_hat_diag_reg',...
        'den_hat', 'den_hat_diag', 'den_hat_diag_reg', 'den_hat_DL', 'den_hat_diag_DL', 'den_hat_diag_reg_DL');
    
    [nr toc]
end
disp('compute_num_sims Finished!');