% Script to run our NRSfM on the Hulk sequence
clear;
close all;
load sceneshulk_additionalviews; % obtained from Chhatkuli et al. data.

totframes = length(scene.m); % number of images
X = [];
Xgth = [];
for k = 1: totframes     
    X = [X;scene.m(k).m(1:2,:)];
    Xgth = [Xgth;scene.Pgth(k).P];
end

idx = isnan(sum(Xgth));
Xt = Xgth(:,~idx);

xmin = min(vec(Xt(1:3:end,:)));
xmax = max(vec(Xt(1:3:end,:)));
ymin = min(vec(Xt(2:3:end,:)));
ymax = max(vec(Xt(2:3:end,:)));
zmin = min(vec(Xt(3:3:end,:)));
zmax = max(vec(Xt(3:3:end,:)));

for k = 1: totframes
    m(k).m = X((2*k-1):(2*k),:);
    Pgt(k).P = Xgth((3*k-2):(3*k),:); % ground truth 3D
end

N = length(m(1).m);
M = length(m);

% visibility is true for the example:
visibt = true(N,M);
Kneighbors = 20;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
IDX = getNeighborsVis(m,Kneighbors,visibt);
visbc = num2cell(visibt,1);

C = getAngleCos(m,IDX);

lambda1 = 1;
lambda2 = 20;

% second part: formulate SDP and solve with cvx
mc  = squeeze(struct2cell(m));
disp('NRSfM function');
tic;
[mu,D] = OurNRSfM(IDX,C,m,lambda1,lambda2);
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
    title('Hulk')    
    clf;
    plot3(Q2k_n(1,:),Q2k_n(2,:),Q2k_n(3,:),'b*');
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    hold on;
    plot3(P2(1,l),P2(2,l),P2(3,l),'go');
    hold off;
    pause(0.2);
    res.Q2{k} = Q2k_n;
    res.Pg{k} = P2(:,l);     
 
    scale = norm(P2(:,l),'fro');    
    res.err3dper(k) = norm(Q2k_n - P2(:,l),'fro')/scale*100;    
    res.err3d(k) = sqrt(mean(sum((P2(:,l)-Q2k_n).^2)));
    fprintf('3D rmse =%.2f mm\t',res.err3d(k));
    fprintf('relative 3D error =%.2f %% \n',res.err3dper(k));
    pause(0.2);        
end
meandepth = mean(res.err3d)
meanper = mean(res.err3dper)

