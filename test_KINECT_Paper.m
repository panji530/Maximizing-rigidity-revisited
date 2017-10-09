% Script to run MDH based NRSfM on the KINECT Paper sequence
clear;
close all;
load KinectPaper; % obtained from Varol et al. data.
load camcalib;
sp = 10; % subsample points 10
totframes = length(p); % number of images
totpoints = length(p(1).p);

Xmat = [];
X = [];
for k = 1: totframes    
    % get normalized points from calibration
    m(k).m = inv(KK)* [p(k).p(1:2,1:sp:end); ones(1,length(p(k).p(:,1:sp:end)))];%Pgth(k).P(:,1:sp:end)./(ones(3,1)*Pgth(k).P(3,1:sp:end));%
    m(k).m = m(k).m(1:2,:);    
    X = [X;m(k).m];
    Pgth(k).P = Pgth(k).P(:,1:sp:end);   
    Xmat = [Xmat;Pgth(k).P];
end

xmin = min(vec(Xmat(1:3:end,:)));
xmax = max(vec(Xmat(1:3:end,:)));
ymin = min(vec(Xmat(2:3:end,:)));
ymax = max(vec(Xmat(2:3:end,:)));
zmin = min(vec(Xmat(3:3:end,:)));
zmax = max(vec(Xmat(3:3:end,:)));

sv = 8; % subsample views
Pgt = Pgth(sv:sv:end);
m = m(sv:sv:end);

N = length(m(1).m);
M = length(m);

% visibility is true for the example:
visibt = true(N,M);
Kneighbors = 20;
lambda1 = 1;
lambda2 = 20;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
[IDX, dist_sort] = getNeighborsVis(m,Kneighbors,visibt);
visbc = num2cell(visibt,1);

C = getAngleCos(m,IDX);

% second part: formulate SDP and solve it with CVX
mc  = squeeze(struct2cell(m));
disp('NRSfM function');
tic;
[mu,D,edgs,maxEdgs] = OurNRSfM(IDX,C, m, lambda1, lambda2);
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
    [Q2k_n,scl] = RegisterToGTH(Q2k(:,l),P2(:,l));  
    
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
