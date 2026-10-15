
% ----- Shape Gradients -----

function [Gx, Gy, Gz, detJ] = shape_gradients7(X, Y, Z, xi)
%   X, Y, Z : (number of elements) x 8 corner coordinates of every element
%   xi      : a point [xi1 xi2 xi3] in the reference cube [0,1]^3
%   Gx, Gy, Gz : dN_a/dx, dN_a/dy, dN_a/dz at that point (elements x 8)
%   detJ       : determinant of the Jacobian for every element

% ----- Step 1 -----
% Derivatives in the reference cube
% Node a sits at corner d = (d1,d2,d3) with each d_i = 0 or 1, and
%   N_a = L1*L2*L3,  where L_i = xi_i if d_i = 1,  or  1 - xi_i if d_i = 0.
dN = zeros(3, 8);
for a = 1:8
    d = [mod(a-1,2), mod(floor((a-1)/2),2), floor((a-1)/4)];
    L  = d.*xi + (1-d).*(1-xi);       % the three 1D linear functions
    dL = 2*d - 1;                     % their derivatives (+1 or -1)
    dN(:,a) = [dL(1)*L(2)*L(3);  L(1)*dL(2)*L(3);  L(1)*L(2)*dL(3)];
end

% ----- Step 2 -----
% The Jacobian J(i,j) = d x_i / d xi_j  (one value per element)
J11 = X*dN(1,:)';  J12 = X*dN(2,:)';  J13 = X*dN(3,:)';
J21 = Y*dN(1,:)';  J22 = Y*dN(2,:)';  J23 = Y*dN(3,:)';
J31 = Z*dN(1,:)';  J32 = Z*dN(2,:)';  J33 = Z*dN(3,:)';

% ----- Step 3 -----
% Chain rule  grad_x N = inv(J)' * grad_xi N
% inv(J)' = (cofactor matrix of J) / det(J), written out for a 3x3 matrix
detJ = J11.*(J22.*J33 - J23.*J32) - J12.*(J21.*J33 - J23.*J31) + J13.*(J21.*J32 - J22.*J31);
C11 =  (J22.*J33 - J23.*J32);  C12 = -(J21.*J33 - J23.*J31);  C13 =  (J21.*J32 - J22.*J31);
C21 = -(J12.*J33 - J13.*J32);  C22 =  (J11.*J33 - J13.*J31);  C23 = -(J11.*J32 - J12.*J31);
C31 =  (J12.*J23 - J13.*J22);  C32 = -(J11.*J23 - J13.*J21);  C33 =  (J11.*J22 - J12.*J21);
Gx = ((C11./detJ)*dN(1,:) + (C12./detJ)*dN(2,:) + (C13./detJ)*dN(3,:));
Gy = ((C21./detJ)*dN(1,:) + (C22./detJ)*dN(2,:) + (C23./detJ)*dN(3,:));
Gz = ((C31./detJ)*dN(1,:) + (C32./detJ)*dN(2,:) + (C33./detJ)*dN(3,:));
end
