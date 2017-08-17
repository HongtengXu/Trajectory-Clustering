function [lambda, R, C] = ParaConfig( TT )

R=[];
C=[];
S=[];
for i=1:length(TT)
    
    traj = TT(i).data;
    speed = traj(2:end,:)-traj(1:end-1,:);
    R=[R; traj(:,2)];
    C=[C; traj(:,1)];
    S=[S;(sqrt(speed(:,1).^2+speed(:,2).^2))];
    
end
R=[min(R),max(R)];
C=[min(C),max(C)];
lambda = var(S);
% rmax=mean(S);
% rmin=min(S);
% 
% wr=rmax:(rmin-rmax)/(Iter-1):rmin;