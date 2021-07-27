for nr = 1:length(NR_vec)
    load(strcat(out_folder,'/outs_',string(nr), '_', string(nq_sim), '_', string(block_ID),'.mat'));
    % ToDo: add block_ID loop and average them
end