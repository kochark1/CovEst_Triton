function compute_num_sims_ZF(path, interm_folder,out_folder, block_ID, TT, ch_samples)
    interm_folder=strcat(path, '/', interm_folder);
    out_folder=strcat(path, '/', out_folder);
    nq_set = [1,6];
    ZF_FLAG = true;
    
    
    nq_index = 1;
    if ZF_FLAG
        for nq_sim = nq_set
            numSimulations_ZF;
            
            nq_index = nq_index + 1;
        end
    end
    
    disp('compute_num_sims_ZF Finished!');

end