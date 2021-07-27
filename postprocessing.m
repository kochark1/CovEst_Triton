function postprocessing()
    interm_folder=strcat(path, '/', interm_folder);
    out_folder=strcat(path, '/', out_folder);
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

        nq_index = nq_index + 1;
    end
    
    disp('postprocessing Finished!');
end