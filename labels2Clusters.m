%% Labels to Connected Components

% clc;
% clear;
% close all;

%% Problem Definition


% load data
load dolphins.mat
data = Problem.A;
data = full(data);

N = size(data , 1);

% adjency matrix
for i=1 : N
    adjM{i} = find(data(i,:)==1);
end

% Objective Function
CostFunction = @(x)fitnessFuncSFLA(community,M,x);

nVar = N;              % Number of Unknown Variables
VarSize = [1 nVar];     % Unknown Variables Matrix Size

VarMax = cellfun(@length,adjM);           % Upper Bound of Unknown Variables
VarMin = arrayfun(@length,VarMax);           % Lower Bound of Unknown Variables

%% Initialization

% Initialize Population Members
Position = arrayfun(@(ind)randi([VarMin(ind), VarMax(ind)],1),1:nVar);
x = bestPos;

%%

clusterNum = 1;
clusterID = zeros(1,length(x));
pos = 1;
posFlagVec = zeros(size(x));
clusterID(pos) = clusterNum;

ngbrMat = zeros(length(x));
colVec = 1:length(x);
colVec(x==0) = [];
rVec = x(x~=0);
rowVec = adjM(sub2ind(size(adjM),rVec,colVec));
indVec = sub2ind(size(ngbrMat),rowVec,colVec);
ngbrMat(indVec) = 1;
ngbrMatSym = ngbrMat | ngbrMat.';

bins = conncomp(digraph(ngbrMat),'Type','weak');

%%
% while ~isempty(posVec)
%     if ~isempty(adjM{pos})
%         posN = adjM{pos}(x(pos));
%         if posN == pos || ~ismember(posN,posVec)
%             clusterNum = clusterNum + 1;
%             pos = posVec(1);
%             continue
%         else
%             posVec(1) = [];
%             pos = posN;
%             clusterID(pos) = clusterNum;
%         end
%     else
%         clusterNum = clusterNum + 1;
%         posVec(1) = [];
%         pos = posVec(1);
%     end
% end

while any(~posFlagVec)
    [posFlagVec,clusterID] = labelerFunc(pos,posFlagVec,clusterID,clusterNum,ngbrMatSym);
    if any(~posFlagVec)
        clusterNum = clusterNum + 1;
        pos = find(~posFlagVec,1,'first');
    end
end

%% Results Check with MATLAB ConnectedComponent Function

results = {'false','true'};
disp(['Results: ',results{(sum(clusterID==bins) == length(x)) + 1}])

%% Functions

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
