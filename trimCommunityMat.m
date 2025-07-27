function trimmedMat = trimCommunityMat(commMat,communities)

trimmedMat = commMat;
for i = 1:size(commMat,2)
    connInd = find(commMat(:,i));
    trimmedMat(connInd(communities(connInd) ~= communities(i)),i) = 0;
end
