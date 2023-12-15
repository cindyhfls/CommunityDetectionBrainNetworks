function g = compute_gamma(W)
% implements optimal gamma computation as in Newman 2016
% Community detection in networks:
% Modularity optimization and maximum likelihood are equivalent
% arXiv:1606.02319v

% Jiaxin Cindy Tu 2023.12.08 - modify i=1:5 to actually check for
% convergence? (in the paper it says it converges in fewer than 10
% iterations I guess Nicolini used 5 arbitrarily? Also it does not work for
% negative edges
% Also it assumes that all groups are statistically similar which is pretty
% unlikley

% First guess of gamma
g = 1;
for i=1:5
    Ci = community_louvain(W,g);
    %length(unique(Ci))
    [~,K,~,~,~,~,Bnorm,~] = comm_mat(W,Ci);
    omega_in = mean(diag(Bnorm));
    omega_out = mean(nonzeros(Bnorm-eye(length(Bnorm)).*Bnorm));
    g = (omega_in - omega_out)/(log(omega_in)-log(omega_out))
end
