function hj = generateh(R_sqr_root,M,L,K,N)
rng shuffle;
h = crandn(M, N, L, K);
for l = 1:L
    for k = 1:K
		hj(:, l, k,:) = (R_sqr_root(:,:,l,k)*h(:,:,l,k)); % (M, L, K, N)
    end
end
end