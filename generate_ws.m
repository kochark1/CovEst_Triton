load(strcat(interm_folder,'/RQWfile.mat'));

if ZF_FLAG
    for kk = 1:number_of_users
        Ws_for_each_k(NR_vec, NQ_vec, nq_set, R_sqrt_root, number_of_cells, number_of_users, number_of_antennas, targetCell, kk, pilotSequenceLength, mu, phi, alpha_R, alpha_Q, Rb, TT, Pb, interm_folder);
    end
else
    Ws_for_each_k(NR_vec, NQ_vec, nq_set, R_sqrt_root, number_of_cells, number_of_users, number_of_antennas, targetCell, targetUser, pilotSequenceLength, mu, phi, alpha_R, alpha_Q, Rb, TT, Pb, interm_folder);
end