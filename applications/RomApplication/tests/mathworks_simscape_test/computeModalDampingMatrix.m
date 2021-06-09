function C = computeModalDampingMatrix(dampingRatio,K,M)

% To avoid numerical issues (such as complex eigenvalues with very small
% imaginary parts), make the matrices exactly symmetric.

    K = (K+K')/2;    % Stiffness matrix
    M = (M+M')/2;    % Mass matrix

% Compute the eigen-decomposition associated with the mass and stiffness
% matrices, sorting the eigenvalues in ascending order and permuting
% the corresponding eigenvectors.

    [V,D] = eig(K,M);
    [d,sortIdxs] = sort(diag(D));
    V = V(:,sortIdxs);

% Due to small numerical errors, the six eigenvalues associated with the
% rigid-body modes may not be exactly zero. To avoid numerical issues,
% check that the first six eigenvalues are close enough to zero. Then
% replace them with exact 0 values.

    assert(all(abs(d(1:6))/abs(d(7)) < 1e-9),'Error due to "zero" eigenvalues.');
    d(1:6) = 0;

% Vectors of generalized masses and natural frequencies

    MV = M*V;
    generalizedMasses = diag(V'*MV);
    naturalFrequencies = sqrt(d);

% Compute the modal damping matrix associated with K and M

    C = MV * diag(2*dampingRatio*naturalFrequencies./generalizedMasses) * MV';

end

