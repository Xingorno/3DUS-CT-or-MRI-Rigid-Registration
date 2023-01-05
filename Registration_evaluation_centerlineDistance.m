% objective: compute the distance between the centerlines
% Inputs: 
%   controlPtsGT: N by 3, the centerline sample points of ground truth
%   controlPtsEst: M by 3, the centerline sample points of estimated segmentation

% Ouput:
%   centerlineDistance1: Euclidean distance from the ground truth centerline to the estimated centerline
%   centerlineDistance2: Euclidean distance from the estimated centerline to the ground 
%   overEst: compute how much over segmentation relative to the ground truth (using length ratio as the metric)
%   underEst: compute how much under segmentation relatie to the ground truth (using length ratio as the metric)
%function [centerlineDistance1, centerlineDistance2, overEst, underEst] = ComputeCenterlineDistance(controlPtsGT, controlPtsEst)

clc
clear all
close all
% Step 1: load in the vessel
%%
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\registration\scene1\MR_centerline\";
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\registration\scene1\MR_centerline\";
dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR\Centerline\";
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\MR\centerline\"
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\CT\centerline\"
[controlPtsGT, brachesIndexGT] = ComputeCenterlineControlPointsAddBlankParts(dirPath);

% HV0502
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\registration\scene1\US_centerline\";
% HV0503
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\registration\scene1\US_centerline\";
% HV0504
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\registration\scene1\US_centerline\";

% HV0202
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\registration\scene1\US_centerline\"
% HV0204
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US4\registration\scene1\US_centerline\"
% HV0203
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US3\registration\scene1\US_centerline\"
% HV0206
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US6\registration\scene1\US_centerline\"

% HV01705
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_05\registration\scene1\US_centerline\"
% % HV01703
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_03\registration\scene1\US_centerline\"
% HV01703
dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_02\registration\scene1\US_centerline\"

% PT013
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\US\registration\scene1\centerline\";
% PT012
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\scene1\centerline\"
[controlPtsEst, branchesIndexEst] = ComputeCenterlineControlPointsAddBlankParts(dirPath);

%% Step2: visualize the branches
colorIndex = {'#0072BD', '#D95319', '#EDB120', 	'#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#0000FF', '#00FFFF', '#FF00FF','#FFFF00', 	'#000000', '#FF0000','#00FFFF', '#EDB120','#0072BD',  '#D95319', '#0072BD','#7E2F8E', '#77AC30', '#D95319', '#EDB120', 	'#7E2F8E', '#0072BD', '#D95319', '#EDB120', '#0072BD', '#D95319', '#EDB120', 	'#7E2F8E', '#77AC30', '#4DBEEE'};
uniqueGTColorIndex = {'#FF0000'};
uniqueEstColorIndex = {'#00FF00'};
[uniqueBranchesGT, iaGT, icGT] = unique(brachesIndexGT, 'stable');
[uniqueBranchesEst, iaEst, icEst] = unique(branchesIndexEst, 'stable');

figure (1)
hold on
for i = 1: size(uniqueBranchesGT, 1)
    branchIndexGT = uniqueBranchesGT(i,1);
    if branchIndexGT < size(uniqueBranchesGT, 1)
        controlPts_temp = controlPtsGT(iaGT(branchIndexGT):iaGT(branchIndexGT+1)-1,:);
    else
        controlPts_temp = controlPtsGT(iaGT(branchIndexGT): end, :);
    end
    if branchIndexGT > 30
        matchIndexGT = branchIndexGT -30;
    else
        matchIndexGT = branchIndexGT;
    end 
    scatter3(controlPts_temp(:,1), controlPts_temp(:,2),controlPts_temp(:,3),'filled','MarkerFaceColor',colorIndex{matchIndexGT});
    text(mean(controlPts_temp(:,1)), mean(controlPts_temp(:,2)),mean(controlPts_temp(:,3)), strcat('GT#',num2str(branchIndexGT)));

end
hold off
title('Ground truth')

figure (2)
hold on
for i = 1: size(uniqueBranchesEst, 1)
    branchIndexEst = uniqueBranchesEst(i);
    if branchIndexEst < size(uniqueBranchesEst, 1)
        controlPts_temp = controlPtsEst(iaEst(branchIndexEst):iaEst(branchIndexEst+1)-1,:);
    else
        controlPts_temp = controlPtsEst(iaEst(branchIndexEst): end, :);
    end
     if branchIndexEst > 30
        matchIndexGT = branchIndexEst -30;
    else
        matchIndexGT = branchIndexEst;
    end 
    scatter3(controlPts_temp(:,1), controlPts_temp(:,2),controlPts_temp(:,3),'filled','MarkerFaceColor',colorIndex{matchIndexGT});
    text(mean(controlPts_temp(:,1)), mean(controlPts_temp(:,2)),mean(controlPts_temp(:,3)), strcat('Est#',num2str(branchIndexEst)));

end
hold off
title('Our segmentation')


figure (3)
hold on
for i = 1: size(uniqueBranchesGT, 1)
    branchIndexGT = uniqueBranchesGT(i,1);
    if branchIndexGT < size(uniqueBranchesGT, 1)
        controlPts_temp = controlPtsGT(iaGT(branchIndexGT):iaGT(branchIndexGT+1)-1,:);
    else
        controlPts_temp = controlPtsGT(iaGT(branchIndexGT): end, :);
    end
    if branchIndexGT > 30
        matchIndexGT = branchIndexGT -30;
    else
        matchIndexGT = branchIndexGT;
    end 
    plot3(controlPts_temp(:,1), controlPts_temp(:,2),controlPts_temp(:,3),'.r','MarkerSize',20);
%     text(mean(controlPts_temp(:,1)), mean(controlPts_temp(:,2)),mean(controlPts_temp(:,3)), strcat('GT#',num2str(branchIndexGT)),'FontSize',14, 'FontWeight', 'bold');

end
% hold off
% title('Ground truth')

% figure (2)
% hold on
for i = 1: size(uniqueBranchesEst, 1)
    branchIndexEst = uniqueBranchesEst(i);
    if branchIndexEst < size(uniqueBranchesEst, 1)
        controlPts_temp = controlPtsEst(iaEst(branchIndexEst):iaEst(branchIndexEst+1)-1,:);
    else
        controlPts_temp = controlPtsEst(iaEst(branchIndexEst): end, :);
    end
    plot3(controlPts_temp(:,1), controlPts_temp(:,2),controlPts_temp(:,3),'.b','MarkerSize',20);
%     text(mean(controlPts_temp(:,1)), mean(controlPts_temp(:,2)),mean(controlPts_temp(:,3)), strcat('Est#',num2str(branchIndexEst)), 'FontSize', 10);

end
% grid on
hold off
title('Our segmentation')

%%

    figure (3)
    plot3(controlPtsGT(:,1),controlPtsGT(:,2), controlPtsGT(:,3), '.r','MarkerSize',20)
    hold on
    plot3(controlPtsEst(:,1),controlPtsEst(:,2), controlPtsEst(:,3), '.b','MarkerSize',20)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds before registration')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')
%% Step 3: combine multiple branches for some complex vessel structures
% Reorder the branches [matching branches + unique branches]

%
% % HV 0502
% branchCombinedIndexGT = [1,4,10];
% branchCombinedIndexEst = [7, 8];  
% branchCombinedIndexGT_1 = [15, 17];
% branchCombinedIndexEst_1 = [10, 12]; 
% branchCombinedIndexGT_2 = [11,12];
% branchCombinedIndexEst_2 = []; 

% % HV 0503
% branchCombinedIndexGT = [15 17 37];
% branchCombinedIndexEst = [10 12 14];  
% branchCombinedIndexGT_1 = [1 4 10 11 12];
% branchCombinedIndexEst_1 = [7 8]; 
% branchCombinedIndexGT_2 = [19 20];
% branchCombinedIndexEst_2 = []; 

% % HV 0504
% branchCombinedIndexGT = [1 4];
% branchCombinedIndexEst = [];  
% branchCombinedIndexGT_1 = [36 42];
% branchCombinedIndexEst_1 = []; 
% branchCombinedIndexGT_2 = [];
% branchCombinedIndexEst_2 = []; 
% branchCombinedIndexGT_3 = [];
% branchCombinedIndexEst_3 = []; 
% 
% % HV 0506
% branchCombinedIndexGT = [15 17];
% branchCombinedIndexEst = [3 4];  
% branchCombinedIndexGT_1 = [4 10];
% branchCombinedIndexEst_1 = [9 11]; 
% branchCombinedIndexGT_2 = [41 43];
% branchCombinedIndexEst_2 = [12 14]; 
% branchCombinedIndexGT_3 = [38 39];
% branchCombinedIndexEst_3 = [15 17]; 

% % HV 0202
% branchCombinedIndexGT = [13 15 17 19];
% branchCombinedIndexEst = [1 2];  
% branchCombinedIndexGT_1 = [7 9];
% branchCombinedIndexEst_1 = [6 7]; 
% branchCombinedIndexGT_2 = [];
% branchCombinedIndexEst_2 = [12 14 15 17]; 
% branchCombinedIndexGT_3 = [];
% branchCombinedIndexEst_3 = [20 21]; 

% % HV 0204
% branchCombinedIndexGT = [1 3];
% branchCombinedIndexEst = [10 12];  
% branchCombinedIndexGT_1 = [26 58];
% branchCombinedIndexEst_1 = []; 
% branchCombinedIndexGT_2 = [21 28];
% branchCombinedIndexEst_2 = []; 
% branchCombinedIndexGT_3 = [];
% branchCombinedIndexEst_3 = []; 

% % HV 0203
% branchCombinedIndexGT = [13 15 17 19];
% branchCombinedIndexEst = [12 14];  
% branchCombinedIndexGT_1 = [6 7 9];
% branchCombinedIndexEst_1 = [21 22 28]; 
% branchCombinedIndexGT_2 = [22 28];
% branchCombinedIndexEst_2 = [4 5]; 
% branchCombinedIndexGT_3 = [];
% branchCombinedIndexEst_3 = [6 7]; 

% % HV 0206
% branchCombinedIndexGT = [11 13 15 17 19];
% branchCombinedIndexEst = [2 4];  
% branchCombinedIndexGT_1 = [28 42 44];
% branchCombinedIndexEst_1 = [6 7]; 
% branchCombinedIndexGT_2 = [25 26];
% branchCombinedIndexEst_2 = [15 17]; 
% branchCombinedIndexGT_3 = [3 6 7];
% branchCombinedIndexEst_3 = []; 


% % HV 017-05
% branchCombinedIndexGT = [19 20 21 22 36];
% branchCombinedIndexEst = [2 3 4 5 11];  
% branchCombinedIndexGT_1 = [29 31 32 33];
% branchCombinedIndexEst_1 = [7 9 10]; 
% branchCombinedIndexGT_2 = [24 26];
% branchCombinedIndexEst_2 = [15 16]; 
% branchCombinedIndexGT_3 = [7 8 9];
% branchCombinedIndexEst_3 = [21 22 25]; 
% 
% % HV 017-03
% branchCombinedIndexGT = [1 4 6 7 8 9 11];
% branchCombinedIndexEst = [1 2 3 5 7];  
% branchCombinedIndexGT_1 = [19 20 21 22];
% branchCombinedIndexEst_1 = [18 28 33 35]; 
% branchCombinedIndexGT_2 = [29 31 32 33];
% branchCombinedIndexEst_2 = [19 21 36]; 
% branchCombinedIndexGT_3 = [];
% branchCombinedIndexEst_3 = [22 24]; 

% % HV 017-02
% branchCombinedIndexGT = [4 6];
% branchCombinedIndexEst = [29 30];  
% branchCombinedIndexGT_1 = [19 36];
% branchCombinedIndexEst_1 = [2 4 10]; 
% branchCombinedIndexGT_2 = [29 31];
% branchCombinedIndexEst_2 = [8 19]; 
% branchCombinedIndexGT_3 = [32 33];
% branchCombinedIndexEst_3 = [22 23]; 
% branchCombinedIndexGT_4 = [];
% branchCombinedIndexEst_4 = [20 25]; 

% % PT013
% branchCombinedIndexGT = [9 10];
% branchCombinedIndexEst = [2 5];  
% branchCombinedIndexGT_1 = [18 21];
% branchCombinedIndexEst_1 = [3 6]; 
% branchCombinedIndexGT_2 = [];
% branchCombinedIndexEst_2 = []; 
% branchCombinedIndexGT_3 = [];
% branchCombinedIndexEst_3 = []; 
% branchCombinedIndexGT_4 = [];
% branchCombinedIndexEst_4 = []; 

% PT012
branchCombinedIndexGT = [37 39];
branchCombinedIndexEst = [18 19];  
branchCombinedIndexGT_1 = [8 9 10 11];
branchCombinedIndexEst_1 = [1 2 3 6]; 
branchCombinedIndexGT_2 = [18 19];
branchCombinedIndexEst_2 = [11 12 13]; 
branchCombinedIndexGT_3 = [49 51 52];
branchCombinedIndexEst_3 = [7 8]; 
branchCombinedIndexGT_4 = [37 39];
branchCombinedIndexEst_4 = [18 19]; 
branchCombinedIndexGT_5 = [33 36 40];
branchCombinedIndexEst_5 = []; 


[controlPtsGT_combined_1, branchesIndexGT_combined_1] = CombineBranches(controlPtsGT, brachesIndexGT, branchCombinedIndexGT);
[controlPtsEst_combined_1, branchesIndexEst_combined_1] = CombineBranches(controlPtsEst, branchesIndexEst, branchCombinedIndexEst);

[controlPtsGT_combined_2, branchesIndexGT_combined_2] = CombineBranches(controlPtsGT_combined_1, branchesIndexGT_combined_1, branchCombinedIndexGT_1);
[controlPtsEst_combined_2, branchesIndexEst_combined_2] = CombineBranches(controlPtsEst_combined_1, branchesIndexEst_combined_1, branchCombinedIndexEst_1);

[controlPtsGT_combined_3, branchesIndexGT_combined_3] = CombineBranches(controlPtsGT_combined_2, branchesIndexGT_combined_2, branchCombinedIndexGT_2);
[controlPtsEst_combined_3, branchesIndexEst_combined_3] = CombineBranches(controlPtsEst_combined_2, branchesIndexEst_combined_2, branchCombinedIndexEst_2);

[controlPtsGT_combined_4, branchesIndexGT_combined_4] = CombineBranches(controlPtsGT_combined_3, branchesIndexGT_combined_3, branchCombinedIndexGT_3);
[controlPtsEst_combined_4, branchesIndexEst_combined_4] = CombineBranches(controlPtsEst_combined_3, branchesIndexEst_combined_3, branchCombinedIndexEst_3);

[controlPtsGT_combined_5, branchesIndexGT_combined_5] = CombineBranches(controlPtsGT_combined_4, branchesIndexGT_combined_4, branchCombinedIndexGT_4);
[controlPtsEst_combined_5, branchesIndexEst_combined_5] = CombineBranches(controlPtsEst_combined_4, branchesIndexEst_combined_4, branchCombinedIndexEst_4);


[controlPtsGT_combined, branchesIndexGT_combined] = CombineBranches(controlPtsGT_combined_5, branchesIndexGT_combined_5, branchCombinedIndexGT_5);
[controlPtsEst_combined, branchesIndexEst_combined] = CombineBranches(controlPtsEst_combined_5, branchesIndexEst_combined_5, branchCombinedIndexEst_5);




% Step 3: re-visualize the vessel branches


% colorIndex = {'#0072BD', '#D95319', '#EDB120', 	'#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', '#0000FF', '#00FFFF', '#FF00FF','#000000','#0072BD', '#D95319', '#EDB120','#0072BD', '#D95319', '#EDB120','#0072BD', '#D95319', '#EDB120' };
uniqueGTColorIndex = {'#FF0000'};
uniqueEstColorIndex = {'#00FF00'};
[uniqueBranchesGT_combined, ~, ~] = unique(branchesIndexGT_combined, 'stable');
[uniqueBranchesEst_combined, ~, ~] = unique(branchesIndexEst_combined, 'stable');
% sortedIaEst = sort(iaEst);
% uniqueBranchesEst_combined_sorted = branchesIndexEst_combined(sortedIaEst);
% sortedIaGT = sort(iaGT);
% uniqueBranchesGT_combined_sorted = branchesIndexGT_combined(sortedIaGT);

% compute the mass center
centerBranchesGT = [];
centerBranchesEst = [];
for i = 1 : size(uniqueBranchesGT_combined, 1)
    branchNumberGT = uniqueBranchesGT_combined(i);
    branchRangeGT = find(branchesIndexGT_combined == branchNumberGT);
    centerBranchesGT = [centerBranchesGT; mean(controlPtsGT_combined(branchRangeGT,:))];
end

for i = 1: size(uniqueBranchesEst_combined, 1)
    branchNumberEst = uniqueBranchesEst_combined(i);
    branchRangeEst = find(branchesIndexEst_combined == branchNumberEst);
    centerBranchesEst = [centerBranchesEst; mean(controlPtsEst_combined(branchRangeEst,:))];
end


% coarse compute the matching braches (from the ground truth to the estiamted segmentation)

branchesMatchingTable = [];
for i = 1: size(centerBranchesGT, 1)
    
    deltaBranches_temp = centerBranchesGT(i, :) - centerBranchesEst;
    deltaBrachesDistance = sqrt(deltaBranches_temp(:,1).*deltaBranches_temp(:,1) + deltaBranches_temp(:,2).*deltaBranches_temp(:,2) + deltaBranches_temp(:,3).*deltaBranches_temp(:,3));
    [minBranchesDistance, minIndex] = min(deltaBrachesDistance);
    if minBranchesDistance < 6
        branchesMatchingTable = [branchesMatchingTable; [uniqueBranchesGT_combined(i), uniqueBranchesEst_combined(minIndex)]]; 
    end
end

%
% branchesMatchingTable(1,:) = [];
% %%%%%%%%%
% HV05_02
% addMatchingTable = [7 11; 11 14; 15 1; 25 4];
% HV05_03
% addMatchingTable = [15, 1; 26, 4; 34 9; 1 10; 7 11; 39 16];
% HV05_04
% addMatchingTable = [1 2; 38 4; 44 5; 45 6; 41 8; 36 9; 15 11];
% HV05_06
% branchesMatchingTable(6,:) = [];
% addMatchingTable = [26 7; 27 8; 41 12; 38 15; 40 16];

% % HV02_02
% branchesMatchingTable(4,:) = [];
% addMatchingTable = [44 3; 45 4; 43 6; 30 11; 13 12; 39 19; 7 20];

% % HV0204
% % branchesMatchingTable(4,:) = [];
% addMatchingTable = [13 4; 26 7; 21 10];
% branchesMatchingTable = [branchesMatchingTable; addMatchingTable];

% % HV0203
% branchesMatchingTable(3:5,:) = [];
% addMatchingTable = [13 12; 3 20; 6 21; 10 25; 23 3; 25 2; 22 4; 50 11];
% branchesMatchingTable = [branchesMatchingTable; addMatchingTable];

% % HV0206
% % branchesMatchingTable(,:) = [];
% addMatchingTable = [12 3; 47 9; 22 11; 25 12; 3 15];
% branchesMatchingTable = [branchesMatchingTable; addMatchingTable];
% %%%%%%%%%%

% % HV17-05

% addMatchingTable = [19 2; 23 6; 35 8; 34 20; 38 19; 39 14; 24 15; 12 17; 16 24; 41 29; 42 30; 43 32];
% branchesMatchingTable = [branchesMatchingTable; addMatchingTable];

% % HV17-03
% branchesMatchingTable(5:6,:) = [];
% addMatchingTable = [1 1; 10 6; 3 8; 13 10; 19 18; 41 13; 42 14];
% branchesMatchingTable = [branchesMatchingTable; addMatchingTable];

% % HV17-02
% branchesMatchingTable(11,:) = [];
% addMatchingTable = [5 46; 32 22; 39 11; 38 14];
% branchesMatchingTable = [branchesMatchingTable; addMatchingTable];

% % PT013
% 
% addMatchingTable = [18 2; 20 7];
% branchesMatchingTable = [branchesMatchingTable; addMatchingTable];

% PT012
branchesMatchingTable(5,:) = [];
addMatchingTable = [13 7; 21 26; 49 24; 33 11; 42 14; 37 18 ];
branchesMatchingTable = [branchesMatchingTable; addMatchingTable];


% branchesMatchingTable

% % branchesMatchingTable(12,2) = 10;
% branchesMatchingTable = [branchesMatchingTable; [13 10]];

% % HV045
% branchesMatchingTable(6,:) = [];
% branchesMatchingTable(4,:) = [];
% 
% % % % HV056 branch
% branchesMatchingTable = [[3 1]; branchesMatchingTable];



branchesNotmatchingGT = [];
branchesNotmatchingEst = [];

sortedMatchingBranchesGT = sort(branchesMatchingTable(:,1),'descend');
branchesNotmatchingGT = uniqueBranchesGT_combined;
for i = 1: size(sortedMatchingBranchesGT, 1)
    branchN = find(branchesNotmatchingGT == sortedMatchingBranchesGT(i));
    branchesNotmatchingGT(branchN) = [];
end


sortedMatchingBranchesEst = sort(branchesMatchingTable(:,2),'descend');
branchesNotmatchingEst = uniqueBranchesEst_combined;
for i = 1 : size(sortedMatchingBranchesEst, 1)
    branchN = find(branchesNotmatchingEst == sortedMatchingBranchesEst(i));
    branchesNotmatchingEst(branchN) = [];
end


figure
hold on
for i = 1: size(branchesMatchingTable, 1)
    branchNumberGT = branchesMatchingTable(i, 1);   
    branchRangeGT = find(branchesIndexGT_combined == branchNumberGT);
    controlPtsGT_temp = controlPtsGT_combined(branchRangeGT,:);

    scatter3(controlPtsGT_temp(:,1), controlPtsGT_temp(:,2),controlPtsGT_temp(:,3),'filled','MarkerFaceColor',colorIndex{i});
    text(mean(controlPtsGT_temp(:,1) + 5), mean(controlPtsGT_temp(:,2)),mean(controlPtsGT_temp(:,3)), strcat('GT#',num2str(branchNumberGT)));
    
    branchNumberEst = branchesMatchingTable(i, 2);
    branchRangeEst = find(branchesIndexEst_combined == branchNumberEst);
    controlPtsEst_temp = controlPtsEst_combined(branchRangeEst,:);

    scatter3(controlPtsEst_temp(:,1), controlPtsEst_temp(:,2),controlPtsEst_temp(:,3),'MarkerEdgeColor',colorIndex{i});
    text(mean(controlPtsEst_temp(:,1)-5), mean(controlPtsEst_temp(:,2)) ,mean(controlPtsEst_temp(:,3)), strcat('Est#', num2str(branchNumberEst)));
end

if ~isempty(branchesNotmatchingGT)
    
    for i = 1:size(branchesNotmatchingGT,1)
        branchNumberGT = branchesNotmatchingGT(i, 1);
        branchRangeGT = find(branchesIndexGT_combined == branchNumberGT);
        controlPtsGT_temp = controlPtsGT_combined(branchRangeGT,:);

        scatter3(controlPtsGT_temp(:,1), controlPtsGT_temp(:,2),controlPtsGT_temp(:,3),'filled','MarkerFaceColor',uniqueGTColorIndex{1});
        text(mean(controlPtsGT_temp(:,1)), mean(controlPtsGT_temp(:,2)),mean(controlPtsGT_temp(:,3)), strcat('GT#',num2str(branchNumberGT)));
    end
end

if ~isempty(branchesNotmatchingEst)
    
    for i = 1:size(branchesNotmatchingEst,1)

         branchNumberEst = branchesNotmatchingEst(i, 1);
         branchRangeEst = find(branchesIndexEst_combined == branchNumberEst);
         controlPtsEst_temp = controlPtsEst_combined(branchRangeEst,:);

         scatter3(controlPtsEst_temp(:,1), controlPtsEst_temp(:,2),controlPtsEst_temp(:,3),'MarkerEdgeColor',uniqueEstColorIndex{1});
         text(mean(controlPtsEst_temp(:,1)), mean(controlPtsEst_temp(:,2)),mean(controlPtsEst_temp(:,3)), strcat('Est#', num2str(branchNumberEst)));
    end
end
hold off

title('Ground truth and estimated segmentation')
legend('GT', 'Est')
%% Step 4: compute the distance

[minDistanceTotalGT2Est,minIndexPairTotalGT2Est, minPointsTotalGT2Est, startingIndexMatchingBranchGT2Est] = ComputeCenterlineDistanceFromGT2Est(controlPtsGT_combined, controlPtsEst_combined, branchesMatchingTable, branchesIndexGT_combined, branchesIndexEst_combined);

controlPtsGT_ = [];
controlPtsEst_ = [];
for k = 1: size(branchesMatchingTable, 1)
    branchNumGT = branchesMatchingTable(k, 1);
    branchNumEst = branchesMatchingTable(k, 2);
    branchIndexSetGT = find(branchesIndexGT_combined == branchNumGT);
    branchIndexSetEst = find(branchesIndexEst_combined == branchNumEst);
    controlPtsGT = controlPtsGT_combined(branchIndexSetGT,:);
    controlPtsEst = controlPtsEst_combined(branchIndexSetEst, :);
    controlPtsGT_ = [controlPtsGT_; controlPtsGT];
    controlPtsEst_ = [controlPtsEst_; controlPtsEst]; 
end

% GT2Est
ptsGT1 = controlPtsGT_combined(minIndexPairTotalGT2Est(:,1),:);
% ptsGT1 = controlPtsGT_(minIndexPairTotalGT2Est(:,1),:);
ptsEst1 = controlPtsEst_combined(minIndexPairTotalGT2Est(:,2),:);
ptsEst1_min = minPointsTotalGT2Est;
%     ptsGT = controlPtsGT_combined(minIndexPairTotal1(:,1),:);
%     ptsEst = controlPtsEst_combined(minIndexPairTotal1(:,2),:);
% xt = [ptsGT1(:,1), ptsEst1_min(:,1)];
% yt = [ptsGT1(:,2), ptsEst1_min(:,2)];
% zt = [ptsGT1(:,3), ptsEst1_min(:,3)];

xt = [ptsGT1(:,1), ptsEst1(:,1)];
yt = [ptsGT1(:,2), ptsEst1(:,2)];
zt = [ptsGT1(:,3), ptsEst1(:,3)];

figure
scatter3(ptsGT1(:,1), ptsGT1(:,2),ptsGT1(:,3),'filled','MarkerFaceColor','#0072BD');
hold on
scatter3(ptsEst1(:,1), ptsEst1(:,2),ptsEst1(:,3), 'MarkerEdgeColor','#D95319');
% scatter3(ptsEst1(:,1),ptsEst1(:,2),ptsEst1(:,3), 'filled', 'r' );

for i = 1:size(ptsGT1,1)
     plot3(xt(i,:), yt(i,:), zt(i,:),'g');
end

if ~isempty(branchesNotmatchingGT)
    for i = 1:size(branchesNotmatchingGT,1)
        branchNumberGT = branchesNotmatchingGT(i, 1);
        branchRangeGT = find(branchesIndexGT_combined == branchNumberGT);
        controlPtsGT_temp = controlPtsGT_combined(branchRangeGT,:);

        scatter3(controlPtsGT_temp(:,1), controlPtsGT_temp(:,2),controlPtsGT_temp(:,3),'filled','MarkerFaceColor',uniqueGTColorIndex{1});
        text(mean(controlPtsGT_temp(:,1)), mean(controlPtsGT_temp(:,2)),mean(controlPtsGT_temp(:,3)), strcat('GT#',num2str(branchNumberGT)));
    end
end

if ~isempty(branchesNotmatchingEst)
    
    for i = 1:size(branchesNotmatchingEst,1)

         branchNumberEst = branchesNotmatchingEst(i, 1);
         branchRangeEst = find(branchesIndexEst_combined == branchNumberEst);
         controlPtsEst_temp = controlPtsEst_combined(branchRangeEst,:);

         scatter3(controlPtsEst_temp(:,1), controlPtsEst_temp(:,2),controlPtsEst_temp(:,3),'MarkerEdgeColor',uniqueEstColorIndex{1});
         text(mean(controlPtsEst_temp(:,1)), mean(controlPtsEst_temp(:,2)),mean(controlPtsEst_temp(:,3)), strcat('Est#', num2str(branchNumberEst)));
    end
end


hold off
legend('Ground truth', 'Estimation')
title("Distance: GT to Est")


[minDistanceTotalEst2GT,minIndexPairTotalEst2GT, minPointsTotalEst2GT, startingIndexMatchingBranchEst2GT] = ComputeCenterlineDistanceFromEst2GT(controlPtsGT_combined, controlPtsEst_combined, branchesMatchingTable, branchesIndexGT_combined, branchesIndexEst_combined);

% Est2GT
figure
ptsGT2 = controlPtsGT_combined(minIndexPairTotalEst2GT(:,1),:);
ptsEst2 = controlPtsEst_combined(minIndexPairTotalEst2GT(:,2),:);
ptsGT2_min = minPointsTotalEst2GT;
%     ptsGT = controlPtsGT_combined(minIndexPairTotal1(:,1),:);
%     ptsEst = controlPtsEst_combined(minIndexPairTotal1(:,2),:);
% xt = [ptsGT2_min(:,1), ptsEst2(:,1)];
% yt = [ptsGT2_min(:,2), ptsEst2(:,2)];
% zt = [ptsGT2_min(:,3), ptsEst2(:,3)];

xt = [ptsGT2(:,1), ptsEst2(:,1)];
yt = [ptsGT2(:,2), ptsEst2(:,2)];
zt = [ptsGT2(:,3), ptsEst2(:,3)];

scatter3(ptsGT2(:,1), ptsGT2(:,2),ptsGT2(:,3),'filled','MarkerFaceColor','#0072BD');
hold on
scatter3(ptsEst2(:,1), ptsEst2(:,2),ptsEst2(:,3), 'MarkerEdgeColor','#D95319');
% scatter3(ptsEst2(:,1),ptsEst2(:,2),ptsEst2(:,3), 'filled', 'r' );

for i = 1:size(ptsGT2,1)
     plot3(xt(i,:), yt(i,:), zt(i,:),'g');
end
    

if ~isempty(branchesNotmatchingGT)
    for i = 1:size(branchesNotmatchingGT,1)
        branchNumberGT = branchesNotmatchingGT(i, 1);
        branchRangeGT = find(branchesIndexGT_combined == branchNumberGT);
        controlPtsGT_temp = controlPtsGT_combined(branchRangeGT,:);

        scatter3(controlPtsGT_temp(:,1), controlPtsGT_temp(:,2),controlPtsGT_temp(:,3),'filled','MarkerFaceColor',uniqueGTColorIndex{1});
        text(mean(controlPtsGT_temp(:,1)), mean(controlPtsGT_temp(:,2)),mean(controlPtsGT_temp(:,3)), strcat('GT#',num2str(branchNumberGT)));
    end
end

if ~isempty(branchesNotmatchingEst)
    
    for i = 1:size(branchesNotmatchingEst,1)

         branchNumberEst = branchesNotmatchingEst(i, 1);
         branchRangeEst = find(branchesIndexEst_combined == branchNumberEst);
         controlPtsEst_temp = controlPtsEst_combined(branchRangeEst,:);

         scatter3(controlPtsEst_temp(:,1), controlPtsEst_temp(:,2),controlPtsEst_temp(:,3),'MarkerEdgeColor',uniqueEstColorIndex{1});
         text(mean(controlPtsEst_temp(:,1)), mean(controlPtsEst_temp(:,2)),mean(controlPtsEst_temp(:,3)), strcat('Est#', num2str(branchNumberEst)));
    end
end

hold off
legend('Ground truth', 'Estimation')
title('Distance: Est to GT')
%% filter those outliers
filterIndexEst2GTTotal = [];
filterIndexGT2EstTotal = []; % in order to remove some noise,we choose three points at the starting/end part to remove


pointNumberGT = 0;
pointNumberEst = 0;
commonVesselLengthEst = [];
commonVesselLengthGT = [];
vesselOverSegmenationEst = [];
vesselLength_over_starting = 0;
vesselLength_over_ending = 0;
vesselUnderSegmenationEst = [];
vesselLength_under_starting = 0;
vesselLength_under_ending = 0;

k =12
% for k = 1: size(branchesMatchingTable, 1)
% for k = 1: 1
    branchNumGT = branchesMatchingTable(k, 1);
    branchNumEst = branchesMatchingTable(k, 2);
    branchIndexSetGT = find(branchesIndexGT_combined == branchNumGT);
    branchIndexSetEst = find(branchesIndexEst_combined == branchNumEst);
    
    startingSide1 = find(minIndexPairTotalEst2GT(:, 1) == min(branchIndexSetGT));
    startingSide2 = find(minIndexPairTotalEst2GT(:, 1) == min(branchIndexSetGT)+ 1);
    startingSide3 = find(minIndexPairTotalEst2GT(:, 1) == min(branchIndexSetGT)+ 2);
    filterIndexEst2GT = [];
    
    if size(startingSide3, 1) > 1 
        filterIndexEst2GT = [filterIndexEst2GT; [startingSide1; startingSide2; startingSide3(1:end-1,:)]];
    elseif size(startingSide2, 1) > 1
        filterIndexEst2GT = [filterIndexEst2GT; [startingSide1; startingSide2(1:end-1,:)]];
    elseif size(startingSide1, 1) > 1 
        filterIndexEst2GT = [filterIndexEst2GT; startingSide1(1:end-1, 1)];
    end
   
   
    endingSide3 = find(minIndexPairTotalEst2GT(:, 1) == max(branchIndexSetGT)-2);
    endingSide2 = find(minIndexPairTotalEst2GT(:, 1) == max(branchIndexSetGT)-1);
    endingSide1 = find(minIndexPairTotalEst2GT(:, 1) == max(branchIndexSetGT));
    
    if size(endingSide3, 1) > 1
        filterIndexEst2GT = [filterIndexEst2GT; [endingSide3(2:end, :); endingSide2; endingSide1]];
    elseif size(endingSide2, 1) > 1
        filterIndexEst2GT = [filterIndexEst2GT; [endingSide2(2:end, :); endingSide1]];
    elseif size(endingSide1, 1) > 1
        filterIndexEst2GT = [filterIndexEst2GT; endingSide1(2:end, :)];
    end
    
    filterIndexEst2GTTotal = [filterIndexEst2GTTotal; filterIndexEst2GT];
    % save out branches
    
    minIndexPairTotalEst2GT_temp  = minIndexPairTotalEst2GT;
    minIndexPairTotalEst2GT_temp(filterIndexEst2GT, :) = [];
    branchIndex = find((minIndexPairTotalEst2GT_temp(:,1)>= min(branchIndexSetGT)) & (minIndexPairTotalEst2GT_temp(:,1)<= max(branchIndexSetGT)));
    branchGT = controlPtsGT_combined(minIndexPairTotalEst2GT_temp(branchIndex,1),:);
    branchEst = controlPtsEst_combined(minIndexPairTotalEst2GT_temp(branchIndex,2),:);
    
    minDistanceTotalEst2GT_temp = minDistanceTotalEst2GT;
    minDistanceTotalEst2GT_temp(filterIndexEst2GT, :) = [];
    minDistanceEst2GT_singleBranch = minDistanceTotalEst2GT_temp(branchIndex, :);
    
    figure(8) 
    scatter3(branchGT(:,1), branchGT(:,2), branchGT(:,3), 'filled','MarkerFaceColor',colorIndex{k});
    hold on
%     pause(0.5)
    scatter3(branchEst(:,1), branchEst(:,2), branchEst(:,3), 'MarkerEdgeColor',colorIndex{k+1});
%     pause(0.5)
    hold off
    title('Est2GT')
    legend('GT','Est')
%     figure (2)
%     plot(minDistanceEst2GT_singleBranch)
%     text(50, 0.5, strcat('Est2GT = ',num2str(mean(minDistanceEst2GT_singleBranch)), '\pm', num2str(std(minDistanceEst2GT_singleBranch))));
%     title('minDistanceEst2GT')
    
    save(strcat('branch',num2str(k), 'Est2GT_GT'), 'branchGT');
    save(strcat('branch',num2str(k), 'Est2GT_Est'), 'branchEst');
    save(strcat('branch',num2str(k), 'Est2GT_minD'), 'minDistanceEst2GT_singleBranch');
%     save('gt2_GT2Est', 'minDistanceTotalGT2Est_filtered');
% end



% filterIndexGT2EstTotal = []; % in order to remove some noise,we choose three points at the starting/end part to remove
% for k = 1: size(branchesMatchingTable, 1)
% for k = 1: 1
    branchNumGT = branchesMatchingTable(k, 1);
    branchNumEst = branchesMatchingTable(k, 2);
    branchIndexSetGT = find(branchesIndexGT_combined == branchNumGT);
    branchIndexSetEst = find(branchesIndexEst_combined == branchNumEst);
%     startingSide = find(minIndexPairTotalGT2Est(:, 2) == min(branchIndexSetEst));
    
    startingSide1 = find(minIndexPairTotalGT2Est(:, 2) == min(branchIndexSetEst));
    startingSide2 = find(minIndexPairTotalGT2Est(:, 2) == min(branchIndexSetEst)+ 1);
    startingSide3 = find(minIndexPairTotalGT2Est(:, 2) == min(branchIndexSetEst)+ 2);
    filterIndexGT2Est = [];
    
    if size(startingSide3, 1) > 1 
        filterIndexGT2Est = [filterIndexGT2Est; [startingSide1; startingSide2; startingSide3(1:end-1,:)]];
    elseif size(startingSide2, 1) > 1
        filterIndexGT2Est = [filterIndexGT2Est; [startingSide1; startingSide2(1:end-1,:)]];
    elseif size(startingSide1, 1) > 1 
        filterIndexGT2Est = [filterIndexGT2Est; startingSide1(1:end-1, 1)];
    end
   
   
    endingSide3 = find(minIndexPairTotalGT2Est(:, 2) == max(branchIndexSetEst)-2);
    endingSide2 = find(minIndexPairTotalGT2Est(:, 2) == max(branchIndexSetEst)-1);
    endingSide1 = find(minIndexPairTotalGT2Est(:, 2) == max(branchIndexSetEst));
    
    if size(endingSide3, 1) > 1
        filterIndexGT2Est = [filterIndexGT2Est; [endingSide3(2:end, :); endingSide2; endingSide1]];
    elseif size(endingSide2, 1) > 1
        filterIndexGT2Est = [filterIndexGT2Est; [endingSide2(2:end, :); endingSide1]];
    elseif size(endingSide1, 1) > 1
        filterIndexGT2Est = [filterIndexGT2Est; endingSide1(2:end, :)];
    end
    
    filterIndexGT2EstTotal = [filterIndexGT2EstTotal; filterIndexGT2Est];
  
    % find the inside part
%     [uniqueElements, ia_, ic_] = unique(minIndexPairTotalGT2Est(branchIndexSetGT,2), 'stable');
    branchIndex_temp = find((minIndexPairTotalGT2Est(:,1)>= min(branchIndexSetGT)) & (minIndexPairTotalGT2Est(:,1)<= max(branchIndexSetGT)));
    for i = (min(branchIndexSetEst)+3) : (max(branchIndexSetEst)-3)
        element_temp = find(minIndexPairTotalGT2Est(:,2) == i)
        if size(element_temp, 1) > 1
            filterIndexGT2Est = [filterIndexGT2Est; element_temp];
        end
    end
    
    
    
    % save out each branch 
    minIndexPairTotalGT2Est_temp  = minIndexPairTotalGT2Est;
    minIndexPairTotalGT2Est_temp(filterIndexGT2Est, :) = [];
    branchIndex = find((minIndexPairTotalGT2Est_temp(:,1)>= min(branchIndexSetGT)) & (minIndexPairTotalGT2Est_temp(:,1)<= max(branchIndexSetGT)));
    branchGT = controlPtsGT_combined(minIndexPairTotalGT2Est_temp(branchIndex,1),:);
    branchEst = controlPtsEst_combined(minIndexPairTotalGT2Est_temp(branchIndex,2),:);
    
    minDistanceTotalGT2Est_temp = minDistanceTotalGT2Est;
    minDistanceTotalGT2Est_temp(filterIndexGT2Est, :) = [];
    minDistanceGT2Est_singleBranch = minDistanceTotalGT2Est_temp(branchIndex, :);
    
    
    
    figure(9)
    
    scatter3(branchGT(:,1), branchGT(:,2), branchGT(:,3), 'filled','MarkerFaceColor',colorIndex{k});
%     pause(0.5)
    hold on
    scatter3(branchEst(:,1), branchEst(:,2), branchEst(:,3), 'MarkerEdgeColor',colorIndex{k+1});
%     pause(0.5)
    hold off
    title('GT2Est')
    legend('GT', 'Est')
    
%     
%     figure (4)
%    
%     plot(minDistanceGT2Est_singleBranch)
%     text(50, 0.5, strcat('GT2Est = ',num2str(mean(minDistanceGT2Est_singleBranch)), '\pm', num2str(std(minDistanceGT2Est_singleBranch))));
%     title('minDistanceGT2Est')

    rangeY = linspace(min(minDistanceGT2Est_singleBranch),  max(minDistanceGT2Est_singleBranch), 5);
%     minDistanceGT2Est_singleBranch(find(minDistanceGT2Est_singleBranch > 3),:) = [];
    figure(10) 
    plot(minDistanceGT2Est_singleBranch)
    text(size(minDistanceGT2Est_singleBranch,1)/2, rangeY(4), strcat('GT2Est = ',num2str(mean(minDistanceGT2Est_singleBranch)), '\pm', num2str(std(minDistanceGT2Est_singleBranch))));
    hold on
    plot(minDistanceEst2GT_singleBranch)
    text(size(minDistanceEst2GT_singleBranch,1)/2, rangeY(2), strcat('Est2GT = ',num2str(mean(minDistanceEst2GT_singleBranch)), '\pm', num2str(std(minDistanceEst2GT_singleBranch))));
    text(size(minDistanceEst2GT_singleBranch,1)/2-10, rangeY(3), strcat('Bidirectional = ',num2str(mean([minDistanceEst2GT_singleBranch; minDistanceGT2Est_singleBranch])), '\pm', num2str(std([minDistanceEst2GT_singleBranch; minDistanceGT2Est_singleBranch]))));
    
    title('minDistanceEst')
    legend('GT2Est', 'Est2GT')
    hold off
    
    save(strcat('branch',num2str(k), 'GT2Est_GT'), 'branchGT');
    save(strcat('branch',num2str(k), 'GT2Est_Est'), 'branchEst');
    save(strcat('branch',num2str(k), 'GT2Est_minD'), 'minDistanceGT2Est_singleBranch');
    
% end


% Step5: compute over segmenation and under segmentation

% pointNumberGT = 0;
% pointNumberEst = 0;
% commonVesselLengthEst = [];
% commonVesselLengthGT = [];
% vesselOverSegmenationEst = [];
% vesselLength_over_starting = 0;
% vesselLength_over_ending = 0;
% vesselUnderSegmenationEst = [];
% vesselLength_under_starting = 0;
% vesselLength_under_ending = 0;

% for k = 1: size(branchesMatchingTable, 1)
    branchNumGT = branchesMatchingTable(k, 1);
    branchNumEst = branchesMatchingTable(k, 2);
    branchIndexSetGT = find(branchesIndexGT_combined == branchNumGT);
    branchIndexSetEst = find(branchesIndexEst_combined == branchNumEst);
    controlPtsGT = controlPtsGT_combined(branchIndexSetGT,:);
    controlPtsEst = controlPtsEst_combined(branchIndexSetEst, :);
    
    vesselLength = 0;
    vesselLengthGT = 0;
    % compute the vessel branch length
    controlPtsEst_shift = controlPtsEst;
    controlPtsEst_shift(1,:) = [];
    controlPtsEst_shift = [controlPtsEst_shift; controlPtsEst(end, :)];
    delta = controlPtsEst_shift - controlPtsEst;
    distance = sqrt(delta(:,1).*delta(:,1) + delta(:, 2).*delta(:,2) + delta(:,3).*delta(:,3));
    noisyPoint = find(distance > 5);
    distance(noisyPoint) = [];
    vesselLength = sum(distance);
    commonVesselLengthEst = [commonVesselLengthEst, vesselLength];
    
 
    % compute the vessel branch length
    controlPtsGT_shift = controlPtsGT;
    controlPtsGT_shift(1,:) = [];
    controlPtsGT_shift = [controlPtsGT_shift; controlPtsGT(end, :)];
    delta = controlPtsGT_shift - controlPtsGT;
    distance = sqrt(delta(:,1).*delta(:,1) + delta(:, 2).*delta(:,2) + delta(:,3).*delta(:,3));
    vesselLengthGT = sum(distance);
    commonVesselLengthGT = [commonVesselLengthGT, vesselLengthGT];
    
%     startingSide = find(minIndexPairTotalEst2GT(:, 1) == min(branchIndexSetGT));
    
    startingSide1 = find(minIndexPairTotalEst2GT(:, 1) == min(branchIndexSetGT));
    startingSide2 = find(minIndexPairTotalEst2GT(:, 1) == min(branchIndexSetGT)+ 1);
    startingSide3 = find(minIndexPairTotalEst2GT(:, 1) == min(branchIndexSetGT)+ 2);
    filterIndexEst2GT_starting = [];
    
    if size(startingSide3, 1) > 1 
        filterIndexEst2GT_starting = [filterIndexEst2GT_starting; [startingSide1; startingSide2; startingSide3(1:end-1,:)]];
    elseif size(startingSide2, 1) > 1
        filterIndexEst2GT_starting = [filterIndexEst2GT_starting; [startingSide1; startingSide2(1:end-1,:)]];
    elseif size(startingSide1, 1) > 1 
        filterIndexEst2GT_starting = [filterIndexEst2GT_starting; startingSide1(1:end-1, 1)];
    end
    
    controlPts_temp = controlPtsEst_combined(minIndexPairTotalEst2GT(filterIndexEst2GT_starting,2), :);
    controlPts_temp_shift = controlPtsEst_combined(minIndexPairTotalEst2GT(filterIndexEst2GT_starting,2), :);
    if isempty(controlPts_temp_shift)
        vesselLength_over_starting= 0;
    else
        controlPts_temp_shift(1,:) = [];
        controlPts_temp_shift = [controlPts_temp_shift; controlPts_temp(end,:)];
        delta = controlPts_temp_shift - controlPts_temp;
        distance = sqrt(delta(:,1).*delta(:,1) + delta(:, 2).*delta(:,2) + delta(:,3).*delta(:,3));
        vesselLength_over_starting = sum(distance);
    end

    

    filterIndexEst2GT_ending = [];
    endingSide3 = find(minIndexPairTotalEst2GT(:, 1) == max(branchIndexSetGT)-2);
    endingSide2 = find(minIndexPairTotalEst2GT(:, 1) == max(branchIndexSetGT)-1);
    endingSide1 = find(minIndexPairTotalEst2GT(:, 1) == max(branchIndexSetGT));
    
    if size(endingSide3, 1) > 1
        filterIndexEst2GT_ending = [filterIndexEst2GT_ending; [endingSide3(2:end, :); endingSide2; endingSide1]];
    elseif size(endingSide2, 1) > 1
        filterIndexEst2GT_ending = [filterIndexEst2GT_ending; [endingSide2(2:end, :); endingSide1]];
    elseif size(endingSide1, 1) > 1
        filterIndexEst2GT_ending = [filterIndexEst2GT_ending; endingSide1(2:end, :)];
    end
    
    controlPts_temp = controlPtsEst_combined(minIndexPairTotalEst2GT(filterIndexEst2GT_ending,2), :);
    controlPts_temp_shift = controlPtsEst_combined(minIndexPairTotalEst2GT(filterIndexEst2GT_ending,2), :);
    if isempty(controlPts_temp_shift)
        vesselLength_over_ending = 0;
    else
        controlPts_temp_shift(1,:) = [];
        controlPts_temp_shift = [controlPts_temp_shift; controlPts_temp(end,:)];
        delta = controlPts_temp_shift - controlPts_temp;
        distance = sqrt(delta(:,1).*delta(:,1) + delta(:, 2).*delta(:,2) + delta(:,3).*delta(:,3));
        vesselLength_over_ending = sum(distance);
    end
    
    
    vesselOverSegmenationEst = [vesselOverSegmenationEst; [vesselLength_over_starting vesselLength_over_ending]];
    
%     controlPtsGT_ = [controlPtsGT_; controlPtsGT];
%     controlPtsEst_ = [controlPtsEst_; controlPtsEst]; 
    disp(strcat('Branch', num2str(branchNumGT)))
    disp('OVER SEGMENTATION: ');
    disp(strcat('starting part: ', num2str(vesselLength_over_starting),'  ending part: ', num2str(vesselLength_over_ending), '  vessel length (Est): ', num2str(vesselLength)));
    disp(strcat('over segmentation ratio: ', num2str((vesselLength_over_starting + vesselLength_over_ending)/vesselLengthGT)));
    
    
    
%     vesselLength = 0;
%     % compute the vessel branch length
%     controlPtsGT_shift = controlPtsGT;
%     controlPtsGT_shift(1,:) = [];
%     controlPtsGT_shift = [controlPtsGT_shift; controlPtsGT(end, :)];
%     delta = controlPtsGT_shift - controlPtsGT;
%     distance = sqrt(delta(:,1).*delta(:,1) + delta(:, 2).*delta(:,2) + delta(:,3).*delta(:,3));
%     vesselLength = sum(distance);
%     commonVesselLengthGT = [commonVesselLengthGT, vesselLength];
     
    
    startingSide1 = find(minIndexPairTotalGT2Est(:, 2) == min(branchIndexSetEst));
    startingSide2 = find(minIndexPairTotalGT2Est(:, 2) == min(branchIndexSetEst)+ 1);
    startingSide3 = find(minIndexPairTotalGT2Est(:, 2) == min(branchIndexSetEst)+ 2);
    filterIndexGT2Est_starting = [];
    
    if size(startingSide3, 1) > 1 
        filterIndexGT2Est_starting = [filterIndexGT2Est_starting; [startingSide1; startingSide2; startingSide3(1:end-1,:)]];
    elseif size(startingSide2, 1) > 1
        filterIndexGT2Est_starting = [filterIndexGT2Est_starting; [startingSide1; startingSide2(1:end-1,:)]];
    elseif size(startingSide1, 1) > 1 
        filterIndexGT2Est_starting = [filterIndexGT2Est_starting; startingSide1(1:end-1, 1)];
    end
    
    controlPts_temp = controlPtsGT_combined(minIndexPairTotalGT2Est(filterIndexGT2Est_starting, 1), :);
    controlPts_temp_shift = controlPtsGT_combined(minIndexPairTotalGT2Est(filterIndexGT2Est_starting, 1), :);
    if isempty(controlPts_temp_shift)
        vesselLength_under_starting = 0;
    else
        controlPts_temp_shift(1,:) = [];
        controlPts_temp_shift = [controlPts_temp_shift; controlPts_temp(end,:)];
        delta = controlPts_temp_shift - controlPts_temp;
        distance = sqrt(delta(:,1).*delta(:,1) + delta(:, 2).*delta(:,2) + delta(:,3).*delta(:,3));
        vesselLength_under_starting = sum(distance);
    end
    
    
    endingSide3 = find(minIndexPairTotalGT2Est(:, 2) == max(branchIndexSetEst)-2);
    endingSide2 = find(minIndexPairTotalGT2Est(:, 2) == max(branchIndexSetEst)-1);
    endingSide1 = find(minIndexPairTotalGT2Est(:, 2) == max(branchIndexSetEst));
    filterIndexGT2Est_ending = [];
    if size(endingSide3, 1) > 1
        filterIndexGT2Est_ending = [filterIndexGT2Est_ending; [endingSide3(2:end, :); endingSide2; endingSide1]];
    elseif size(endingSide2, 1) > 1
        filterIndexGT2Est_ending = [filterIndexGT2Est_ending; [endingSide2(2:end, :); endingSide1]];
    elseif size(endingSide1, 1) > 1
        filterIndexGT2Est_ending = [filterIndexGT2Est_ending; endingSide1(2:end, :)];
    end
    
    
    controlPts_temp = controlPtsGT_combined(minIndexPairTotalGT2Est(filterIndexGT2Est_ending, 1), :);
    controlPts_temp_shift = controlPtsGT_combined(minIndexPairTotalGT2Est(filterIndexGT2Est_ending, 1), :);
    if isempty(controlPts_temp_shift)
        vesselLength_under_ending = 0;
    else
        controlPts_temp_shift(1,:) = [];
        controlPts_temp_shift = [controlPts_temp_shift; controlPts_temp(end,:)];
        delta = controlPts_temp_shift - controlPts_temp;
        distance = sqrt(delta(:,1).*delta(:,1) + delta(:, 2).*delta(:,2) + delta(:,3).*delta(:,3));
        vesselLength_under_ending = sum(distance);
    end

    vesselUnderSegmenationEst = [vesselUnderSegmenationEst; [vesselLength_under_starting vesselLength_under_ending]];
    
%     controlPtsGT_ = [controlPtsGT_; controlPtsGT];
%     controlPtsEst_ = [controlPtsEst_; controlPtsEst]; 
%     disp(strcat('UNDER SEGMENTATION: ','Branch', num2str(branchNumGT)));
    disp('UNDER SEGMENTATION: ');
    disp(strcat('starting part: ', num2str(vesselLength_under_starting),'  ending part: ', num2str(vesselLength_under_ending), '  vessel length (GT): ', num2str(vesselLengthGT)));
    disp(strcat('under segmentation ratio: ', num2str((vesselLength_under_starting + vesselLength_under_ending)/vesselLengthGT)));
    
    
    % display original vessel
    figure(11)
    scatter3(controlPtsGT(:,1), controlPtsGT(:,2), controlPtsGT(:,3), 'filled','MarkerFaceColor',colorIndex{k});
    hold on
    scatter3(controlPtsEst(:,1), controlPtsEst(:,2), controlPtsEst(:,3), 'MarkerEdgeColor',colorIndex{k+1});
    hold off
    legend('GT', 'Est')
    title('orginal vessel branch')
% end

%%
% compute filtered matching vessel index

minDistanceTotalGT2Est_filtered = minDistanceTotalGT2Est;
minDistanceTotalGT2Est_filtered(filterIndexGT2EstTotal) = [];

minDistanceTotalEst2GT_filtered = minDistanceTotalEst2GT;
minDistanceTotalEst2GT_filtered(filterIndexEst2GTTotal) = [];


figure (6)
plot(minDistanceTotalGT2Est_filtered)
text(100, 0.15, strcat('GT2Est = ',num2str(mean(minDistanceTotalGT2Est_filtered)), '\pm', num2str(std(minDistanceTotalGT2Est_filtered))));
hold on
% mean(minDistanceTotalGT2Est_filtered)
plot(minDistanceTotalEst2GT_filtered)
text(50, 0.1, strcat('Est2GT = ',num2str(mean(minDistanceTotalEst2GT_filtered)), '\pm', num2str(std(minDistanceTotalEst2GT_filtered))));
% mean(minDistanceTotalEst2GT_filtered)
text(200, 0.4, strcat('Total = ',num2str(mean([minDistanceTotalGT2Est_filtered; minDistanceTotalEst2GT_filtered])), '\pm', num2str(std([minDistanceTotalGT2Est_filtered; minDistanceTotalEst2GT_filtered]))));
hold off
legend('GT2Est', 'Est2GT')

%%
save('gt2_GT2Est', 'minDistanceTotalGT2Est_filtered');
save('gt2_Est2GT', 'minDistanceTotalEst2GT_filtered');





