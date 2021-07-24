function [R_est, Q_est]= getRQest(R_sqr_root, L, K, M, targetCell, targetUser, pilotSequenceLength, NR, NQ,  mu, phi, alpha_R, Rb)
N1 = crandn(M, pilotSequenceLength, NR);
N2 = crandn(M, pilotSequenceLength, NR);
phiConj = conj(phi(:, targetCell, targetUser)); % conj(p_{u}), conjugate of the target user ChEst pilot.

R_est = 0;
hj = generateh(R_sqr_root,M,L,K,NR); % generate channel vectors for all the users to the BS j (targetCell), here hj is M X NR X L X K tensor.
%hj(:, n, l, k) is the channel h_{jlk}[n] in the paper.
for n = 1:NR
    epslion1 = 0;
    for l = 1:L
        if l ~= targetCell
            epslion1 = epslion1 + hj(:,l,targetUser,n); % pilot contmaniation in eq 7 
        end
    end
%     epslion2 = epslion1*exp(-1i*2*pi*rand);
    epslion2 = epslion1*exp(1i*2*pi*(rand-rand)); % pilot contmaniation in eq 8. here both rands corresponds to theta_{ln} and theta_{jn}.
    % Addtional note: the assumption on knowledge of theta_{jn} is already
    % present in obtaining eq 9 so no need to explicitly handle it in the
    % code.

    h1 = hj(:,targetCell, targetUser,n) + epslion1 + (1/(pilotSequenceLength*sqrt(mu)))*N1(:,:,n)*phiConj; % eq 7 in paper
%     h2 = hj(:,targetCell, targetUser,n) + epslion2 + (1/(pilotSequenceLength*sqrt(mu)))*N2(:,:,n)*phiConj;
    h2 = hj(:,targetCell, targetUser,n) + epslion2 + (1/(pilotSequenceLength*sqrt(mu)))*N2(:,:,n)*phiConj*exp(-1i*2*pi*rand); % eq 8 in paper
    % Here, phiConj*exp(-1i*2*pi*rand) is the  conjugate of the target user CovEst pilot.
    
    
%     h1 = hj(:,targetCell, targetUser,n); % eq 7 in paper
%     h2 = hj(:,targetCell, targetUser,n); % eq 7 in paper

    R_est = R_est + (1/(2*NR))*(h1*h2' + h2*h1');
end
R_est = alpha_R*R_est+(1-alpha_R)*Rb;

% Similarly, we generate Q matrix from a different set of NQ ChEst pilots.
Q_est = 0;
N = crandn(M, pilotSequenceLength, NQ);
hj = generateh(R_sqr_root,M,L,K,NQ);
for n = 1:NQ
    hsum = 0;
    for l = 1:L
        hsum = hsum + hj(:,l,targetUser,n);
    end
    h_hat = hsum + (1/(pilotSequenceLength*sqrt(mu)))*N(:,:,n)*phiConj;    
    Q_est = Q_est + (1/NQ)*(h_hat*h_hat');
end

end