% ----- Assemble K -----

function K = assemble_K8(M, S, sigma_enc)
%   K = assemble_K(M, S, sigma_enc): tissue uses the 3x3 tensor S and the
%       encapsulation elements use sigma_enc * I.
%   K = assemble_K(M, S): S everywhere (used in the verification tests).
%
%   For each element:  K_e(a,b) = integral of grad(N_a)' * sigma * grad(N_b) dV
%   using 2x2x2 Gauss-Legendre quadrature.  All elements are done at once.

useEnc = (nargin == 3);
X = M.X(M.conn);  Y = M.Y(M.conn);  Z = M.Z(M.conn);    
nE = size(M.conn, 1);
Ke = zeros(nE, 64);           
gp = [0.5 - 0.5/sqrt(3), 0.5 + 0.5/sqrt(3)];   % Gauss points on [0,1]
w = 1/8;                                        % each of the 8 points has weight 1/8 on [0,1]^3

for p = 1:2
  for q = 1:2
    for s = 1:2
        [Gx, Gy, Gz, detJ] = shape_gradients7(X, Y, Z, [gp(p) gp(q) gp(s)]);
        % sigma * grad N
        SGx = S(1,1)*Gx + S(1,2)*Gy + S(1,3)*Gz;
        SGy = S(2,1)*Gx + S(2,2)*Gy + S(2,3)*Gz;
        SGz = S(3,1)*Gx + S(3,2)*Gy + S(3,3)*Gz;
        if useEnc                                % encapsulation: scalar conductivity
            SGx(M.enc,:) = sigma_enc*Gx(M.enc,:);
            SGy(M.enc,:) = sigma_enc*Gy(M.enc,:);
            SGz(M.enc,:) = sigma_enc*Gz(M.enc,:);
        end
        wdet = w*abs(detJ);                      % quadrature weight * volume factor
        for b = 1:8
            for a = 1:8
                col = (b-1)*8 + a;
                Ke(:,col) = Ke(:,col) + wdet .* (Gx(:,a).*SGx(:,b) + Gy(:,a).*SGy(:,b) + Gz(:,a).*SGz(:,b));
            end
        end
    end
  end
end

% Add every element matrix into the big sparse matrix (sparse() adds repeated entries)
rows = M.conn(:, repmat(1:8, 1, 8));          % node of a   for column (b-1)*8 + a
cols = M.conn(:, kron(1:8, ones(1,8)));       % node of b
N = length(M.X);
K = sparse(rows(:), cols(:), Ke(:), N, N);
K = (K + K')/2;                                % exactly symmetric
end
