
dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\US\registration\centerlineDistance\"
fileExtension_min = "*minD.mat";
minCenterlineDistance = [];

% for j = 1:13
%     branchname = strcat(dirPath, 'gt', num2str(j), '\');
    branchname = dirPath;
    searchingFiles = strcat(branchname, fileExtension_min);
    matFiles = dir(searchingFiles);
    

for i = 1: size(matFiles, 1)
%     i = 1;
    nameTemp = matFiles(i).name;
    fullName = strcat(branchname, nameTemp)
    load(fullName);
    
    %a.minDistanceEst2GT_singleBranch;
    
    minCenterlineDistance = [minCenterlineDistance; minDistanceEst2GT_singleBranch;];
    
end
% end
figure 
plot(minCenterlineDistance)

mean(minCenterlineDistance)
std(minCenterlineDistance)
minCenterlineDistance10 = minCenterlineDistance;

%%
fileExtension_pts_Est = '*Est2GT_Est.mat'
fileExtension_pts_GT = '*Est2GT_GT.mat'
branchname = dirPath;
searchingFiles_pts_Est = strcat(branchname, fileExtension_pts_Est);
searchingFiles_pts_GT = strcat(branchname, fileExtension_pts_GT);
matFiles_pts_Est = dir(searchingFiles_pts_Est);
matFiles_pts_GT = dir(searchingFiles_pts_GT);
figure
hold
for i = 1: size(matFiles_pts_Est, 1)
%     i = 1;
    nameTemp_Est = matFiles_pts_Est(i).name;
    
    fullName_Est = strcat(branchname, nameTemp_Est);
    
    nameTemp_GT = matFiles_pts_GT(i).name;
    fullName_GT = strcat(branchname, nameTemp_GT)
    
    load(fullName_Est);
    load(fullName_GT);
    
    %a.minDistanceEst2GT_singleBranch;
    scatter3(branchEst(:,1), branchEst(:,2), branchEst(:,3), 'r')
    scatter3(branchGT(:,1), branchGT(:,2), branchGT(:,3), 'g')
%     minCenterlineDistance = [minCenterlineDistance; minDistanceEst2GT_singleBranch;];
    
end
