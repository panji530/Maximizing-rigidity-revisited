% Script to run our NRSfM on the T-Shirt sequence
clear;
close all;
load tshirtdata; % obtained from Chhatkuli et al. data.

totframes = length(scene.m); % number of images
X = [];
Xgth = [];
for k = 1: totframes        
    X = [X;scene.m(k).m(1:2,:)];
    Xgth = [Xgth;scene.Pgth(k).P];
end

idx = isnan(sum(Xgth));
X = X(:,~idx);
Xgth = Xgth(:,~idx);

xmin = min(vec(Xgth(1:3:end,:)));
xmax = max(vec(Xgth(1:3:end,:)));
ymin = min(vec(Xgth(2:3:end,:)));
ymax = max(vec(Xgth(2:3:end,:)));
zmin = min(vec(Xgth(3:3:end,:)));
zmax = max(vec(Xgth(3:3:end,:)));

for k = 1: totframes
    m(k).m = X((2*k-1):(2*k),:);
    Pgt(k).P = Xgth((3*k-2):(3*k),:);
end

N = length(m(1).m);
M = length(m);

% visibility is true for the example:
visibt = true(N,M);
Kneighbors = 20;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
[IDX, dist_sort] = getNeighborsVis(m,Kneighbors,visibt);
visbc = num2cell(visibt,1);

C = getAngleCos(m,IDX);

lambda1 = 1;
lambda2 = 20;

% second part: formulate SDP and solve with CVX
mc  = squeeze(struct2cell(m));
disp('NRSfM function');
tic;
[mu,D] = OurNRSfM(IDX,C, m, lambda1,lambda2);
ts = toc;

% third part: display results, 
%%
res.Q2 = cell(1,M);
res.Pg = cell(1,M);
res.err3d = zeros(1,M); % RMSE for each surface
res.err3dper = zeros(1,M); % RMSE for each surface
for k=1:M
    Q2k=double([mu(k,visibt(:,k));mu(k,visibt(:,k));mu(k,visibt(:,k))]).*[m(k).m(:,visibt(:,k));ones(1,length(m(k).m(:,visibt(:,k))))];    
    P2 = Pgt(k).P(:,visibt(:,k));
    % get valid indices: some groundtruth points are 0
    mugth = P2(3,:);
    l = mugth>0;
    % fix scale of reconstructed surface
    Q2k_n = RegisterToGTH(Q2k(:,l),P2(:,l));
    
    figure(1)
    clf;    
    plot3(Q2k_n(1,:),Q2k_n(2,:),Q2k_n(3,:),'b*');
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    hold on;
    plot3(P2(1,l),P2(2,l),P2(3,l),'go');
    hold off;
    pause(0.1);
    res.Q2{k} = Q2k_n;
    res.Pg{k} = P2(:,l);   
    
 
    scale = norm(P2(:,l),'fro');    
    res.err3dper(k) = norm(Q2k_n - P2(:,l),'fro')/scale*100;    
    res.err3d(k) = sqrt(mean(sum((P2(:,l)-Q2k_n).^2)));
    fprintf('3D rmse =%.2f mm\t',res.err3d(k));
    fprintf('relative 3D error =%.2f %% \n',res.err3dper(k));       
end
meandepth = mean(res.err3d)
meanper = mean(res.err3dper)

