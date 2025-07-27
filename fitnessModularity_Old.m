function [ Energy ] = fitnessModularity(x, community ,adjM)
%FITNESS_FUNC Summary of this function goes here
%   Detailed explanation goes here


colVec = 1:length(x);
colVec(x==0) = [];
rowVec = x(x~=0);
L=zeros(size(x));
L(x~=0) = adjM(sub2ind(size(adjM),rowVec,colVec));

n_community = length(unique(community));
for i =1 : n_community
    per_community{i} = find(community==i);
end

M = zeros(length(L));
indVec = sub2ind(size(M),L(x~=0),colVec);
M(indVec) = 1;
 
% for k = 1 : length(L)
%     if L(1,k) > 0
%         M(k,L(1,k)) = 1;
%     end
% end

total_edge = sum(sum(M));

% Q = 0;
% for k = 1 : n_community
%     % if ~isempty(per_community{k})
%         edge_community = length(per_community{k}) - ( length(per_community{k}) - length(intersect( per_community{k}',L(1,per_community{k}) ) ));
%         degree_community = sum(sum( M(per_community{k},:) ));
%         Q = Q + ( (edge_community / total_edge) - (degree_community / (2*total_edge)).^2 );
%     % end
% end
% Energy = -Q;

edge_community = cellfun(@(comm)length(comm) - (length(comm) - length(intersect(comm',L(1,comm)))),per_community);
degree_community = cellfun(@(comm)sum(sum(M(comm,:))),per_community);
Q = arrayfun(@(x,y)(x/total_edge) - (y/(2*total_edge)).^2,edge_community,degree_community);

Energy = -sum(Q);

end

