% This work was done while the first author was in the University of Adelaide. Copyright reserved.
function C = getAngleCos(m, NgIdx)
% compute the negative cosine of angles between legs
% Inputs: m -- 2D image coordinates
%         NgIdx -- index of neighbors of each point
%         KK -- intrinsic camera matrix
% Output: C -- matrix containing the consine of angles
% Author: Pan Ji, University of Adelaide
F = length(m); % number of frames
N = length(m(1).m); % number of points
%Sn = size(NgIdx,2);
C = struct([]);

for f = 1:F
    tmp_C = eye(N,N);
    X = m(f).m; % normalised image coordinates of the current frame
    X_bar = [X;ones(1,N)]; % add an all-one row
    for i = 1:N
        NgX = X_bar(:,NgIdx(i,:));
        NormOfNgX = sqrt(sum(NgX.^2));
        Xi = X_bar(:,i);
        tmp_C(i,NgIdx(i,:)) = -(Xi'*NgX)./(norm(Xi).*NormOfNgX);
    end
    tmp_C = makeSymmetricMat(tmp_C);
    tmp_C = tmp_C-diag(diag(tmp_C))+eye(N);
    C(f).c = sparse(tmp_C);
end

end
