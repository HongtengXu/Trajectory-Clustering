function [tildeTT,IND] = FastAMKS( TT, wr, lambda, RR, CC, Map )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adaptive Multi-Kernel-based trajectory Shrinkage (AMKS)
%
% Hongteng Xu
% Oct. 4, 2015
% Georgia Tech
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmin=floor(min(RR)); rmax=ceil(max(RR));
cmin=floor(min(CC)); cmax=ceil(max(CC));



IND=zeros(rmax-rmin+1,cmax-cmin+1);

NumTraj = length(TT);
NumPoint = 0;
num=zeros(NumTraj+1,1);
for i=1:NumTraj
    num(i+1) = size(TT(i).data,1);
    NumPoint=NumPoint + num(i+1);
end

Points = zeros(2, NumPoint);
Speed = zeros(size(Points));
Start = zeros(size(Points));
Stop = zeros(size(Points));


for i=1:NumTraj
    traj = TT(i).data;
    for j=1:size(traj,1)
        R = max([rmin, floor(traj(j,2))]);
        C = max([cmin, floor(traj(j,1))]);
        IND(R-rmin+1,C-cmin+1)=IND(R-rmin+1,C-cmin+1)+1;
    end
    xtmp=traj(2:end,:)-traj(1:end-1,:);
    speed=[xtmp;xtmp(end,:)];
    start=traj(1,:);
    stop=traj(end,:);

    lb = 1+sum(num(1:i));
    ub = sum(num(1:i+1));
    Points(:,lb:ub)=traj';
    Speed(:,lb:ub)=speed';
    Start(:,lb:ub)=repmat(start',[1,num(i+1)]);
    Stop(:,lb:ub)=repmat(stop',[1,num(i+1)]);
end

G=fspecial('gaussian', [ceil(wr)*2+1,ceil(wr)*2+1], wr/3);
IND=conv2(log(IND+1),G,'same');

hatTT = zeros(size(Speed));
hatSpeed = zeros(size(Speed));

% tic
for i=1:NumPoint
    points = Points(:,Map{i});
    speeds = Speed(:,Map{i});
    start = Start(:,Map{i});
    stop = Stop(:,Map{i});
    res = sqrt( (Points(1,i)-points(1,:)).^2+(Points(2,i)-points(2,:)).^2 );

    dis1 = res;
    R=floor(points(2,:));
    R(R<rmin)=rmin;
    C=floor(points(1,:));
    C(C<cmin)=cmin;
    dis4=zeros(1,length(R));
    for j=1:length(R);
        dis4(j)=IND(R(j)-rmin+1,C(j)-cmin+1);
    end
    
    
    dis3 = sqrt( (Speed(1,i)-speeds(1,:)).^2 + (Speed(2,i)-speeds(2,:)).^2 );
    dis2 = sqrt( (Start(1,i)-start(1,:)).^2 +...
                (Start(2,i)-start(2,:)).^2 +...
                (Stop(1,i)-stop(1,:)).^2 +...
                (Stop(2,i)-stop(2,:)).^2 );
            
    Weight = (dis4.^1).*exp(-(dis3.^2)./(wr^2)).*exp(-(dis1.^2)./(wr^2)).*exp(-dis2./wr);
    Weight = Weight./sum(Weight);
    
    hatTT(:,i) = points*Weight(:);
    hatSpeed(:,i) = speeds*Weight(:);
    

    if mod(i,10000)==0
        fprintf('FastSIS: %d/%d Points, time=%0.2f sec\r', i, NumPoint, toc);
    end
end

tildeTT=TT;

% min_X ||X-Y||_F^2 + \lambda*||W(LX-Z)||_F^2

% tic;
% options=optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');%,...
for i=1:NumTraj
    lb = 1+sum(num(1:i));
    ub = sum(num(1:i+1));
    now = hatTT(:,lb:ub)';
    speed_now = now(2:end,:)-now(1:end-1,:);
    
    Time = num(i+1);
    L=-[eye(Time-1),zeros(Time-1,1)]+[zeros(Time-1,1),eye(Time-1)];
    I=eye(Time);
    
    W1 = diag(sum(abs(speed_now).^2,2));
    W2 = I;%diag(max(sum(abs(now-ori),2))-sum(abs(now-ori),2));
    
    
    
    A = W2'*W2 + lambda*L'*(W1'*W1)*L;
    
    Y = now;
    Z = hatSpeed(:,lb:ub-1)';%speed_ori; 
    B = lambda*L'*(W1'*W1)*Z+(W2'*W2)*Y;
    X = A\B;
    
    tildeTT(i).data=X;
    

    if mod(i,1000)==0
        fprintf('FastAMKS: %d/%d Trajectories, time=%0.2f sec\r', i, NumTraj, toc);
    end
end
