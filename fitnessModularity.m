function [ Energy ] = fitnessModularity(x,adjM)

community = adjacencyLableToClusterID(x,adjM);

M = zeros(length(x));
colVec = repmat(1:length(x),size(adjM,1),1);
indVec = sub2ind(size(M),adjM(adjM~=0),colVec(adjM~=0));
M(indVec) = 1;

deltaMat = (community == community.');
totalEdges = sum(M(:));
Q = 1/totalEdges*(M - sum(M,1).*(sum(M,1).')/totalEdges).*deltaMat;

Energy = -sum(Q(:));

end

