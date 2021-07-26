function compute_num_sims(path, interm_folder,out_folder, block_ID, TT)
    interm_folder=strcat(path, '/', interm_folder);
    out_folder=strcat(path, '/', out_folder);
    nq_set = [1,6];
    
    nq_index = 1;
    for nq_sim = nq_set
        numSimulations;

        nq_index = nq_index + 1;
    end
    
    disp('compute_num_sims Finished!');

end