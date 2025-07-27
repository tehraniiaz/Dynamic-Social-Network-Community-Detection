function clusterID = adjacencyLableToClusterID(x,adjM)

clusterNum = 1;
clusterID = zeros(size(x));
pos = 1;
posFlagVec = zeros(size(x));
clusterID(pos) = clusterNum;

ngbrMat = zeros(length(x));
colVec = 1:length(x);
colVec(x==0) = [];
rVec = x(x~=0);
rowVec = adjM(sub2ind(size(adjM),rVec,colVec));
% rowVec = arrayfun(@(k)adjM{k}(x(k)),colVec);
indVec = sub2ind(size(ngbrMat),rowVec(rowVec~=0),colVec(rowVec~=0));
ngbrMat(indVec) = 1;
ngbrMat = ngbrMat | ngbrMat.';

while any(~posFlagVec)
    [posFlagVec,clusterID] = labelerFunc(pos,posFlagVec,clusterID,clusterNum,ngbrMat);
    if any(~posFlagVec)
        clusterNum = clusterNum + 1;
        pos = find(~posFlagVec,1,'first');
    end
end

function [posFlagVec,clusterID] = labelerFunc(pos,posFlagVec,clusterID,clusterNum,ngbrMat)

clusterID(pos) = clusterNum;
posN = find(ngbrMat(:,pos));

if isempty(posN)
    posFlagVec(pos) = 1;
    return
end

for i = 1:length(posN)
    if ~posFlagVec(posN(i))
        clusterID(posN(i)) = clusterNum;
        posFlagVec(posN(i)) = 1;
        [posFlagVec,clusterID] = labelerFunc(posN(i),posFlagVec,clusterID,clusterNum,ngbrMat);
    else
        continue
    end
end

end
end