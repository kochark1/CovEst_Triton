number_of_block = 10;
for nr = 1:length(NR_vec)
    
    num_hat_avg = 0;
    num_hat_diag_avg = 0;
    num_hat_diag_reg_avg = 0;

    den_hat_avg = 0;
    den_hat_diag_avg = 0;
    den_hat_diag_reg_avg = 0;
    
    
    den_hat_DL_avg = 0;
    den_hat_diag_DL_avg = 0;
    den_hat_diag_reg_DL_avg = 0;
    
    for block_ID = 1:number_of_block
        load(strcat(out_folder,'/outs_',string(nr), '_', string(nq_sim), '_', string(block_ID),'.mat'));
        num_hat_avg = num_hat_avg + (1/number_of_block)*num_hat;
        num_hat_diag_avg = num_hat_diag_avg + (1/number_of_block)*num_hat_diag;
        num_hat_diag_reg_avg = num_hat_diag_reg_avg + (1/number_of_block)*num_hat_diag_reg;

        den_hat_avg = den_hat_avg + (1/number_of_block)*den_hat;
        den_hat_diag_avg = den_hat_diag_avg + (1/number_of_block)*den_hat_diag;
        den_hat_diag_reg_avg = den_hat_diag_reg_avg + (1/number_of_block)*den_hat_diag_reg;
        
        den_hat_DL_avg = den_hat_DL_avg + (1/number_of_block)*den_hat_DL;
        den_hat_diag_DL_avg = den_hat_diag_DL_avg + (1/number_of_block)*den_hat_diag_DL;
        den_hat_diag_reg_DL_avg = den_hat_diag_reg_DL_avg + (1/number_of_block)*den_hat_diag_reg_DL;
    end

    gamma_hat = real(num_hat_avg/den_hat_avg);
    gamma_hat_diag = real(num_hat_diag_avg/den_hat_diag_avg);
    gamma_hat_diag_reg = real(num_hat_diag_reg_avg/den_hat_diag_reg_avg);

    SE_hat(nr) = prelogFactor(nr)*log2(1+gamma_hat);
    SE_hat_diag(nr) = prelogFactor(nr)*log2(1+gamma_hat_diag);
    SE_hat_diag_reg(nr) = prelogFactor(nr)*log2(1+gamma_hat_diag_reg);

    gamma_hat_DL = real(num_hat_avg/den_hat_DL_avg);
    gamma_hat_diag_DL = real(num_hat_diag_avg/den_hat_diag_DL_avg);
    gamma_hat_diag_reg_DL = real(num_hat_diag_reg_avg/den_hat_diag_reg_DL_avg);

    SE_hat_DL(nr) = log2(1+gamma_hat_DL);
    SE_hat_diag_DL(nr) = log2(1+gamma_hat_diag_DL);
    SE_hat_diag_reg_DL(nr) = log2(1+gamma_hat_diag_reg_DL);
end