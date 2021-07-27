function postprocessing(path, interm_folder,out_folder)
    interm_folder=strcat(path, '/', interm_folder);
    out_folder=strcat(path, '/', out_folder);
    results_folder = 'resultsSE/';
    if ~exist(strcat(path, '/', results_folder), 'dir')
       mkdir(path, results_folder)
    end
    
    nq_set = [1,6];
    ZF_FLAG = true;
    load(strcat(interm_folder,'/RQWfile.mat'));
    
    % Do Not move the following above the param_save function call
    ulResults = UlResults(nq_set);
    ulResults_diag = UlResults(nq_set);
    ulResults_diag_reg = UlResults(nq_set);

    dlResults = DlResults(nq_set);
    dlResults_diag = DlResults(nq_set);
    dlResults_diag_reg = DlResults(nq_set);

    if ZF_FLAG
        ulResults_ZF = UlResults(nq_set);
        ulResults_diag_ZF = UlResults(nq_set);
        ulResults_diag_reg_ZF = UlResults(nq_set);

        dlResults_ZF = DlResults(nq_set);
        dlResults_diag_ZF = DlResults(nq_set);
        dlResults_diag_reg_ZF = DlResults(nq_set);
    end

    load_true_and_theo;

    prelogFactor_LMMSE = (1 - pilotSequenceLength/ulCoherenceSymbols);

    nq_index = 1;
    for nq_sim = nq_set
        load_num_sim;
        if ZF_FLAG
            load_num_sim_ZF;
        end
        
        ulResults(nq_index).trueSE = prelogFactor_LMMSE*SE_const*ones(size(SE_theo(:,nq_sim)));
        ulResults_diag(nq_index).trueSE = prelogFactor_LMMSE*SE_const_diag*ones(size(SE_theo_diag(:,nq_sim)));
        ulResults_diag_reg(nq_index).trueSE = prelogFactor_LMMSE*SE_const_diag*ones(size(SE_theo_diag_reg(:,nq_sim)));

        ulResults(nq_index).theoreticalSE = SE_theo(:,nq_sim);
        ulResults_diag(nq_index).theoreticalSE = SE_theo_diag(:,nq_sim);
        ulResults_diag_reg(nq_index).theoreticalSE = SE_theo_diag_reg(:,nq_sim);

        ulResults(nq_index).numericalSE = SE_hat;
        ulResults_diag(nq_index).numericalSE = SE_hat_diag;
        ulResults_diag_reg(nq_index).numericalSE = SE_hat_diag_reg;

        %-------------
        dlResults(nq_index).trueSE = SE_const_DL*ones(size(SE_theo_DL(:,nq_sim)));
        dlResults_diag(nq_index).trueSE = SE_const_diag_DL*ones(size(SE_theo_diag_DL(:,nq_sim)));
        dlResults_diag_reg(nq_index).trueSE = SE_const_diag_DL*ones(size(SE_theo_diag_reg_DL(:,nq_sim)));

        dlResults(nq_index).theoreticalSE = SE_theo_DL(:,nq_sim);
        dlResults_diag(nq_index).theoreticalSE = SE_theo_diag_DL(:,nq_sim);
        dlResults_diag_reg(nq_index).theoreticalSE = SE_theo_diag_reg_DL(:,nq_sim);

        dlResults(nq_index).numericalSE = SE_hat_DL;
        dlResults_diag(nq_index).numericalSE = SE_hat_diag_DL;
        dlResults_diag_reg(nq_index).numericalSE = SE_hat_diag_reg_DL;

        if ZF_FLAG
            ulResults_ZF(nq_index).numericalSE = SE_hat_ZF;
            ulResults_diag_ZF(nq_index).numericalSE = SE_hat_diag_ZF;
            ulResults_diag_reg_ZF(nq_index).numericalSE = SE_hat_diag_reg_ZF;

            dlResults_ZF(nq_index).numericalSE = SE_hat_DL_ZF;
            dlResults_diag_ZF(nq_index).numericalSE = SE_hat_diag_DL_ZF;
            dlResults_diag_reg_ZF(nq_index).numericalSE = SE_hat_diag_reg_DL_ZF;
        end

        nq_index = nq_index + 1;
    end
    
    save(strcat(results_folder,'uldlResults.mat'), 'NR_vec', 'ulResults', 'ulResults_diag', 'ulResults_diag_reg', 'dlResults', 'dlResults_diag', 'dlResults_diag_reg', 'N_thr_UL','N_thr_DL', 'ulResults_ZF', 'ulResults_diag_ZF', 'ulResults_diag_reg_ZF', 'dlResults_ZF', 'dlResults_diag_ZF', 'dlResults_diag_reg_ZF', 'ZF_FLAG');

    disp('postprocessing Finished!');
end
