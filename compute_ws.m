function compute_ws(path, interm_folder,out_folder, k_input, TT)
    interm_folder=strcat(path, '/', interm_folder);
    out_folder=strcat(path, '/', out_folder);
    load(strcat(interm_folder,'/RQWfile.mat'));
    Ws_for_each_k(NR_vec, NQ_vec, nq_set, R_sqrt_root, number_of_cells,...
        number_of_users, number_of_antennas, targetCell, k_input,...
        pilotSequenceLength, mu_val, phi, alpha_R, alpha_Q, Rb, TT, Pb,...
        interm_folder);
end