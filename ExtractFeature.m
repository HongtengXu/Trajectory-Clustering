function [Traj, Feature, truth] = ExtractFeature(TT, D, flag1, flag2, lr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extract feature via Smoothing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
truth = zeros(length(TT),1);

G=fspecial('gaussian',[2*lr+1,1],(lr/3));

if flag1 == 0 % without re-sampling
    Traj=TT;
    L = size(TT(1).data,1)*size(TT(1).data,2);
    Feature = zeros(L,length(TT));
    for i=1:length(TT)
        traj = TT(i).data;
        if flag2==1
            traj = conv2(traj,G,'valid');
        end
        if isempty(TT(i).label)
            truth(i)=0;
        else
            truth(i) = TT(i).label;
        end
        Feature(1:length(traj(:)),i)=traj(:);
    end
    
else
    Traj=TT;
    
    Feature = zeros(2*D,length(TT));
    for i=1:length(TT)
        traj = TT(i).data;
        
        if isempty(TT(i).label)
            truth(i)=0;
        else
            truth(i) = TT(i).label;
        end
        x=imresize(traj(:,1),[D,1]);
        y=imresize(traj(:,2),[D,1]);
        if flag2==1
            x = conv2(x,G,'valid');
            y = conv2(y,G,'valid');
        end
        Feature(1:length(x)*2,i)=[x;y];
        Traj(i).data=[x,y];
    end
end

