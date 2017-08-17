%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comparison for various trajectory clustering methods, include:
% 1. mean shift (MS)
% 2. manifold blurring mean shift (MBMS) in 
%   "Manifold blurring mean shift algorithms for manifold denoising",
%   Wang, Weiran and Carreira-Perpin{\'a}n, Miguel A,
%   CVPR, 2010
% 3. Our AMKS without speed regularization
% 4. Our AMKS
%
% This code is only for reseach study. Please cite our paper:
%   Unsupervised Trajectory Clustering via Adaptive Multi-Kernel-based Shrinkage,
%   Hongteng Xu, Yang Zhou, Weiyao Lin, Hongyuan Zha
%   ICCV 2015
%
% Hongteng Xu
% Oct. 4, 2015
% Georgia Tech
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

load TrajData.mat
load color.mat




Traj=DataHu2;
Traj1 = DataHu2;
Traj2 = DataHu2;
Traj3 = DataHu2;
Traj4 = DataHu2;

Iter = 20;
rmax = 12;
rmin = 1;
wr=rmax:(rmin-rmax)/(Iter-1):rmin;
lambda=0.5;
lr=5;

flag=1;
flagf=0;
gr=3;
D=50;
[Traj1, ~, ~] = ExtractFeature(Traj1, D, flag, flagf, gr);
[Traj2, ~, ~] = ExtractFeature(Traj2, D, flag, flagf, gr);
[Traj3, ~, ~] = ExtractFeature(Traj3, D, flag, flagf, gr);
[Traj4, ~, ~] = ExtractFeature(Traj4, D, flag, flagf, gr);

figure
subplot(3,5,11)
img1=imread('DataHu2.png');
image(img1)
hold on
for i=1:length(Traj)
    traj=Traj(i).data;
    label=Traj(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight

for k=1:Iter
    if k==1  
        [Traj1, Map] = MeanShift( Traj1, wr(k) );
        Traj2 = MBMSFast( Traj2, wr(k), lr, Map );
        [~, R, C] = ParaConfig( Traj3 );
        [Traj3,~] = FastAMKS( Traj3, wr(k), 0, R, C, Map );
        [~, R, C] = ParaConfig( Traj4 );
        [Traj4,~] = FastAMKS( Traj4, wr(k), lambda, R, C, Map );
    else
        Traj1 = MeanShiftFast( Traj1, wr(k), Map );
        Traj2 = MBMSFast( Traj2, wr(k), lr, Map );
        [~, R, C] = ParaConfig( Traj3 );
        [Traj3,~] = FastAMKS( Traj3, wr(k), 0, R, C, Map );
        [~, R, C] = ParaConfig( Traj4 );
        [Traj4,~] = FastAMKS( Traj4, wr(k), lambda, R, C, Map );
    end
end


subplot(3,5,12)
image(img1)
hold on
for i=1:length(Traj1)
    traj=Traj1(i).data;
    label=Traj1(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('MeanShift');

subplot(3,5,13)
image(img1)
hold on
for i=1:length(Traj2)
    traj=Traj2(i).data;
    label=Traj2(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('MBMS');

subplot(3,5,14)
image(img1)
hold on
for i=1:length(Traj3)
    traj=Traj3(i).data;
    label=Traj3(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('Without Speed Regularization')

subplot(3,5,15)
image(img1)
hold on
for i=1:length(Traj4)
    traj=Traj4(i).data;
    label=Traj4(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('Proposed FastAMKS')


%%
Traj=DataLin;
Traj1 = DataLin;
Traj2 = DataLin;
Traj3 = DataLin;
Traj4 = DataLin;

Iter = 7;
rmax = 45;
rmin = 10;
wr=rmax:(rmin-rmax)/(Iter-1):rmin;
lambda=0.05;%*ones(Iter,1);
lr=10;


flag=1;
flagf=0;
gr=3;
D=50;
[Traj1, ~, ~] = ExtractFeature(Traj1, D, flag, flagf, gr);
[Traj2, ~, ~] = ExtractFeature(Traj2, D, flag, flagf, gr);
[Traj3, ~, ~] = ExtractFeature(Traj3, D, flag, flagf, gr);
[Traj4, ~, ~] = ExtractFeature(Traj4, D, flag, flagf, gr);


subplot(3,5,6)
img2=imread('DataLin.png');
image(img2)
hold on
for i=1:length(Traj)
    traj=Traj(i).data;
    label=Traj(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight

for k=1:Iter
    if k==1  
        [Traj1, Map] = MeanShift( Traj1, wr(k) );
        Traj2 = MBMSFast( Traj2, wr(k), lr, Map );
        [~, R, C] = ParaConfig( Traj3 );
        [Traj3,~] = FastAMKS( Traj3, wr(k), 0, R, C, Map );
        [~, R, C] = ParaConfig( Traj4 );
        [Traj4,~] = FastAMKS( Traj4, wr(k), lambda, R, C, Map );
    else
        Traj1 = MeanShiftFast( Traj1, wr(k), Map );
        Traj2 = MBMSFast( Traj2, wr(k), lr, Map );
        [~, R, C] = ParaConfig( Traj3 );
        [Traj3,~] = FastAMKS( Traj3, wr(k), 0, R, C, Map );
        [~, R, C] = ParaConfig( Traj4 );
        [Traj4,~] = FastAMKS( Traj4, wr(k), lambda, R, C, Map );
    end
end


subplot(3,5,7)
image(img2)
hold on
for i=1:length(Traj1)
    traj=Traj1(i).data;
    label=Traj1(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('MeanShift');

subplot(3,5,8)
image(img2)
hold on
for i=1:length(Traj2)
    traj=Traj2(i).data;
    label=Traj2(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('MBMS');

subplot(3,5,9)
image(img2)
hold on
for i=1:length(Traj3)
    traj=Traj3(i).data;
    label=Traj3(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('Without Speed Regularization')

subplot(3,5,10)
image(img2)
hold on
for i=1:length(Traj4)
    traj=Traj4(i).data;
    label=Traj4(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('Proposed FastAMKS')



%%
Traj=DataMorris;
Traj1 = DataMorris;
Traj2 = DataMorris;
Traj3 = DataMorris;
Traj4 = DataMorris;

Iter = 7;
rmax = 3;
rmin = 1;
wr=rmax:(rmin-rmax)/(Iter-1):rmin;
lambda=10;
lr=10;


flag=1;
flagf=0;
gr=3;
D=50;
[Traj1, ~, ~] = ExtractFeature(Traj1, D, flag, flagf, gr);
[Traj2, ~, ~] = ExtractFeature(Traj2, D, flag, flagf, gr);
[Traj3, ~, ~] = ExtractFeature(Traj3, D, flag, flagf, gr);
[Traj4, ~, ~] = ExtractFeature(Traj4, D, flag, flagf, gr);

subplot(3,5,1)
img3=imread('DataMorris.png');
image(img3)
hold on
for i=1:length(Traj)
    traj=Traj(i).data;
    label=Traj(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight

for k=1:Iter
    if k==1  
        [Traj1, Map] = MeanShift( Traj1, wr(k) );
        Traj2 = MBMSFast( Traj2, wr(k), lr, Map );
        [~, R, C] = ParaConfig( Traj3 );
        [Traj3,~] = FastAMKS( Traj3, wr(k), 0, R, C, Map );
        [~, R, C] = ParaConfig( Traj4 );
        [Traj4,~] = FastAMKS( Traj4, wr(k), lambda, R, C, Map );
    else
        Traj1 = MeanShiftFast( Traj1, wr(k), Map );
        Traj2 = MBMSFast( Traj2, wr(k), lr, Map );
        [~, R, C] = ParaConfig( Traj3 );
        [Traj3,~] = FastAMKS( Traj3, wr(k), 0, R, C, Map );
        [~, R, C] = ParaConfig( Traj4 );
        [Traj4,~] = FastAMKS( Traj4, wr(k), lambda, R, C, Map );
    end
end


subplot(3,5,2)
image(img3)
hold on
for i=1:length(Traj1)
    traj=Traj1(i).data;
    label=Traj1(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('MeanShift');

subplot(3,5,3)
image(img3)
hold on
for i=1:length(Traj2)
    traj=Traj2(i).data;
    label=Traj2(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('MBMS');

subplot(3,5,4)
image(img3)
hold on
for i=1:length(Traj3)
    traj=Traj3(i).data;
    label=Traj3(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('Without Speed Regularization')

subplot(3,5,5)
image(img3)
hold on
for i=1:length(Traj4)
    traj=Traj4(i).data;
    label=Traj4(i).label;
    plot(traj(:,1),traj(:,2),'color',color(label,:));
end
hold off
axis tight
title('Proposed FastAMKS')