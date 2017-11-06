% This work was done while the first author was in the University of Adelaide. Copyright reserved.
function [ IDX,dist_sort ] = getNeighborsVis( m, Ng, visb, opt )
% Author: Ajad Chhatkuli et al.
% Modified by Pan Ji
if(nargin<4)
    opt = 'cityblock';
end
% get triangulation using pdist functions: Ng number of neighbors
N = length(m(1).m);
distmat = zeros(size(m(1).m,2),size(m(1).m,2),length(m));
for k =1: length(m)
    distmat(:,:,k) = pdist2(m(k).m',m(k).m',opt);    
    distmat(~visb(:,k),:,k) = -1;   

end

dist = max(distmat,[],3);

[dist_sort, IDX] = sort(dist,2);

IDX = IDX(:,1:Ng);
dist_sort = dist_sort(:,1:Ng);

end
