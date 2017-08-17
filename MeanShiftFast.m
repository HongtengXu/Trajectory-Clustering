function tildeTT = MeanShiftFast( TT, wr, Map )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Manifold blurring mean shift
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NumTraj = length(TT);
NumPoint = 0;
num=zeros(NumTraj+1,1);
for i=1:NumTraj
    num(i+1) = size(TT(i).data,1);
    NumPoint=NumPoint + num(i+1);
end

Points = zeros(2, NumPoint);



for i=1:NumTraj
    traj = TT(i).data;
    
    lb = 1+sum(num(1:i));
    ub = sum(num(1:i+1));
    Points(:,lb:ub)=traj';
    
end



hatTT = zeros(size(Points));

tic
for i=1:NumPoint
    ind=Map{i};
    points = Points(:,ind);
    res = sqrt( (Points(1,i)-points(1,:)).^2+(Points(2,i)-points(2,:)).^2 );

    
    dis1 = res;
    
    
    Weight = exp(-(dis1.^2)./(wr^2));
    Weight = Weight./sum(Weight);
    
    hatTT(:,i) = points*Weight(:);
    

    if mod(i,10000)==0
        fprintf('MSFast: %d/%d Points, time=%0.2f sec\r', i, NumPoint, toc);
    end
end

tildeTT=TT;

% min_X ||X-Y||_F^2 + \lambda*||W(LX-Z)||_F^2

tic;
% options=optimoptions(@fminunc,'Algorithm','quasi-newton','Display','off');%,...

% LLL=2;
% for nn=1:LLL
for i=1:NumTraj
    lb = 1+sum(num(1:i));
    ub = sum(num(1:i+1));
    now = hatTT(:,lb:ub)';
    
    tildeTT(i).data=now;
   
    if mod(i,1000)==0
        fprintf('MSFast: %d/%d Trajectories, time=%0.2f sec\r', i, NumTraj, toc);
    end

% end
end