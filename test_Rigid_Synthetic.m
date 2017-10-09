clear,
close all
load KinectPaper;

M = 20;
sp = 25;
X = Pgth(100).P(:,1:sp:end);
N = length(X);
Xgth = zeros(3,N,M);
Xgth(:,:,1) = X;
figure(1)
plot3(X(1,:),X(2,:),X(3,:),'bo')

pi = 3.1415926;

Xk = X;
for k=1:M-1    
    theta = 0.5*3.141592*k/(2*M);
    R = [1,0,0;0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];    
    theta = k*pi/(M);
    t = 400*[cos(theta);1;sin(theta)];
    Xk = R*X + t*ones(1,N);
    Xgth(:,:,k+1) = Xk;    
end

xmin = min(vec(Xgth(1:3:end,:)));
xmax = max(vec(Xgth(1:3:end,:)));
ymin = min(vec(Xgth(2:3:end,:)));
ymax = max(vec(Xgth(2:3:end,:)));
zmin = min(vec(Xgth(3:3:end,:)));
zmax = max(vec(Xgth(3:3:end,:)));

for k =1:M
    Xk = Xgth(:,:,k);
    figure(1)
    clf;    
    plot3(Xk(1,:),Xk(2,:),Xk(3,:),'bo')    
    axis equal
    axis([xmin xmax ymin ymax zmin zmax])
    pause(0.1)  
end

for k = 1:M
    m(k).m = Xgth(:,:,k)./(ones(3,1)*Xgth(3,:,k));
    m(k).m = m(k).m(1:2,:);    
    Pgt(k).P = Xgth(:,:,k);     
end

% visibility is true for the example:
visibt = true(N,M);
Kneighbors = 20;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
[IDX, ~] = getNeighborsVis(m,Kneighbors,visibt);
visbc = num2cell(visibt,1);

C = getAngleCos(m,IDX);

% second part: formulate SDP and solve
mc  = squeeze(struct2cell(m));
disp('NRSfM function');
lambda1 = 1;
lambda2 = 20;
tic;
[mu,lgs,edgs,maxEdgs] = OurNRSfM(IDX,C, m, lambda1,lambda2);
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
    l = 1:length(mugth);%mugth>0;
    % fix scale of reconstructed surface
    [Q2k_n, scl] = RegisterToGTH(Q2k(:,l),P2(:,l));
    %scl
    
    figure(2)
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
