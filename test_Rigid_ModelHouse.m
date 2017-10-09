clear,
close all
load model_house;
% script to run our method on the rigid sequence Model_House
% (http://www.robots.ox.ac.uk/~vgg/data/data-mview.html)

xmin = min(vec(Xgt(1:3:end,:)));
xmax = max(vec(Xgt(1:3:end,:)));
ymin = min(vec(Xgt(2:3:end,:)));
ymax = max(vec(Xgt(2:3:end,:)));
zmin = min(vec(Xgt(3:3:end,:)));
zmax = max(vec(Xgt(3:3:end,:)));

[D,N] = size(W);
M = D/2;

for k = 1:M
    m(k).m = inv(KK(k).K)*[W(2*k-1:2*k,:);ones(1,N)];
    m(k).m = m(k).m(1:2,:);
    Pgt(k).P = Xgt;
end

visibt = true(N,M);
Kneighbors = 20;
lambda1 = 1;
lambda2 = 20;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
IDX = getNeighborsVis(m,Kneighbors,visibt);
visbc = num2cell(visibt,1);

C = getAngleCos(m,IDX);

% second part: formulate SDP and solve with CVX
mc  = squeeze(struct2cell(m));
disp('NRSfM function');
tic;
[mu,D,edgs,maxEdgs] = OurNRSfM(IDX,C, m, lambda1, lambda2);
ts = toc;

% third part: display results, 
%%
hlres.Q2 = cell(1,M);
hlres.Pg = cell(1,M);
hlres.err3d = zeros(1,M); % RMSE for each surface
hlres.err3dper = zeros(1,M); % RMSE for each surface
for k=1:M
    Q2k=double([mu(k,visibt(:,k));mu(k,visibt(:,k));mu(k,visibt(:,k))]).*[m(k).m(:,visibt(:,k));ones(1,length(m(k).m(:,visibt(:,k))))];    
    P2 = Pgt(k).P(:,visibt(:,k));
    % get valid indices: some groundtruth points are 0
    mugth = P2(3,:);
    l = 1:N;
    % fix scale of reconstructed surface    
    [~,Q2k_n] = procrustes(P2',Q2k');
    Q2k_n = Q2k_n';
    
    figure(1)
    clf;
    plot3(Q2k_n(1,:),Q2k_n(2,:),Q2k_n(3,:),'b*');
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    hold on;
    plot3(P2(1,l),P2(2,l),P2(3,l),'go');
    hold off;
    pause(0.1);        
    
    hlres.Q2{k} = Q2k_n;
    hlres.Pg{k} = P2(:,l);
 
    scale = norm(P2(:,l),'fro');    
    hlres.err3dper(k) = norm(Q2k_n - P2(:,l),'fro')/scale*100;    
    hlres.err3d(k) = sqrt(mean(sum((P2(:,l)-Q2k_n).^2)));
    fprintf('3D rmse =%.2f mm\t',hlres.err3d(k));
    fprintf('relative 3D error =%.2f %% \n',hlres.err3dper(k));       
end
meandepth = mean(hlres.err3d)
meanper = mean(hlres.err3dper)

