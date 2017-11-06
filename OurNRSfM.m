function [z, lgs, edgs, maxEdgs] = OurNRSfM(IDX, C, m, lambda1, lambda2)
% Author: Pan Ji, University of Adelaide
%         pan.ji@adelaide.edu.au
%         All rights reserved
% Input:  IDX -- K-NNG neighboring point index
%         C   -- matrix of cosine of viewing angles
%         m   -- normalized image measurement
% Output: z   -- recovered depth
%         lgs -- recovered length of legs (see our paper)
%       edges -- recovered edges (see our paper)
%     maxEdgs -- recovered edges in the internal model (maximum distance model)
% Reference:
% Pan Ji, Hongdong Li, Yuchao Dai, Ian Reid. "Maximizing rigidity"
% revisited: a convex programming approach for generic 3D shape
% reconstruction from multiple perspective views. in ICCV 2017.
if(nargin<5)
    lambda2 = 20;
end
if(nargin<4)
    lambda1 = 1;
end

IDX = IDX(:,2:end);
F = length(C); % number of frames
[N, k]= size(IDX); % number of points, number of neighbors

tmpI = speye(N);
A = struct([]);
tmp = sparse(N*k,N*N);

for f = 1:F
    for i = 1:N
        for j = 1:k
            e = diag(tmpI(:,i)+tmpI(:,IDX(i,j)));
            eCe = e*C(f).c'*e;
            eCe = eCe(:);
            tmp((i-1)*k+j,:) = eCe;
        end
    end
    A(f).vecA = tmp;
end

cvx_begin
    cvx_precision low 
    variables x(N*F) d(N*k) l(N*k,F) Y(N,N,F);      
    tmp = trace(Y(:,:,1));
    for i = 2:F
        tmp = tmp+trace(Y(:,:,i));
    end    
    minimize( tmp-lambda1*sum(x)-lambda2*sum(vec(l)) )    
    subject to
    sum(d) == 1;
    x >= 0;
    l >= 0;    
    l <= d*ones(1,F);      
    for f = 1:F    
        idx1 = N*(f-1)+1;
        idx2 = N*f;
        [1, x(idx1:idx2)';x(idx1:idx2), Y(:,:,f)] == semidefinite(N+1);          
        (A(f).vecA)*vec(Y(:,:,f)) == l(:,f);
    end
cvx_end


lgs = x;
edgs = l;
maxEdgs = d;

lgs = reshape(lgs, N, F);
lgs = lgs';
z = zeros(F,N);

for i = 1:F
    qi = [m(i).m; ones(1,N)];
    z(i,:) = lgs(i,:)./sqrt(sum(qi.^2));
end

end
