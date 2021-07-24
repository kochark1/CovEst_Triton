clear all;
number_of_cells = 7; % L
number_of_users = 10; % K
number_of_antennas = 100; % M
ulCoherenceSymbols = 100; %C_u
fullFramecoherenceBlocks = 25000;


targetCell = 1; % j
targetUser = 1; % u

NR_vec = 2.^(-2:3)*1000;
NQ_vec = 2.^(-2:3)*1000;

% NR_vec = [100,500,600];
% NQ_vec = [100,500,600];

mu = 1; % UL signal power
lambda_DL = 10; %DL signal power
pilotSequenceLength = number_of_users; %p
alpha_R = 0.95;
Rb = eye(number_of_antennas);
alpha_Q = 0.95;
Pb = diag(diag(Rb));

% number_of_cells, number_of_users, number_of_antennas, targetCell, targetUser, ULresources, NR, rho, lambda, sigma_sqr, mu

[Rmatrices, Rsum, Rsum_DL, R_sqrt_root, Qmatrix, Wmatrix, W_bar, SE_const, SE_const_DL, W_diag, W_bar_diag, SE_const_diag, SE_const_diag_DL] = ConfigRQWC(number_of_cells, number_of_users, number_of_antennas, targetCell, targetUser, mu, lambda_DL, alpha_R, Rb, pilotSequenceLength);
SE_sum = 0;
SE_sum_diag = 0;
for  targetUser = 1:number_of_users
    [Rmatrices, Rsum, Rsum_DL, R_sqrt_root, Qmatrix, Wmatrix, W_bar, SE_const, SE_const_DL, W_diag, W_bar_diag, SE_const_diag, SE_const_diag_DL] = ConfigRQWC(number_of_cells, number_of_users, number_of_antennas, targetCell, targetUser, mu, lambda_DL, alpha_R, Rb, pilotSequenceLength);
    SE_sum = SE_sum + SE_const;
    SE_sum_diag = SE_sum_diag + SE_const_diag;
    
end
SE_sum
SE_sum_diag
targetUser = 1;