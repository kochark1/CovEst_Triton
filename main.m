close all;
clear all;
clc;

t_start = datetime
tic;
interm_folder = './interm';
out_folder = './out_folder';


TT = 2;
nq_set = [1,6];
ZF_FLAG = true;

% param_save;
% generate_ws;

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

trueAndTheoretical;

prelogFactor_LMMSE = (1 - pilotSequenceLength/ulCoherenceSymbols);

nq_index = 1;
for nq_sim = nq_set
    numSimulations;
    if ZF_FLAG
        numSimulations_ZF;
    end
    
    nq_index = nq_index + 1;
end



toc;
t_end = datetime
between(t_start,t_end)