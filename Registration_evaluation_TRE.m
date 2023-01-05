clc;
clear all;
close all;


%% Step 1: Read vessel centerlines
    % HV-05
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\registration\scene1\MR_centerline\"
    
    % HV-02
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Centerline\";
    
%     % HV017
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR\Centerline\";
%     % PT013
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\MR\centerline\";    
    % PT012
    dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\CT\centerline\";  
    
    fileExtension = "*.json";
    searchingFiles = strcat(dirPath, fileExtension);
    jsonFiles = dir(searchingFiles);
    allCTControlPoints = [];
    allCTBifurcationPoints = [];

    % resample the vessel points for registration
    for i = 1:size(jsonFiles, 1)
    % for i = 8:10
    % for i = 10:size(jsonFiles, 1)
        fileHeader = jsonFiles(i);
        fileName = fileHeader.name;
        fileFolder = fileHeader.folder;
        fileAbsoluteName = strcat(fileFolder, "\", fileName);
        fid = fopen(fileAbsoluteName); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        val = jsondecode(str);
        controlPoints = val.markups.controlPoints;
        allCTControlPoints_temp = [];
        for j = 1:size(controlPoints, 1)
            temp = controlPoints(j);
            allCTControlPoints_temp = [allCTControlPoints_temp; temp.position'];
        end
    %     scaling = ceil(0.2*exp(allCTRadius(i)));% Add weights to stablize registration
        scaling = 1;
        numOfInter = scaling*size(allCTControlPoints_temp,1);
        % interpolate using parametric splines
        pt = interparc(numOfInter,allCTControlPoints_temp(:,1),allCTControlPoints_temp(:,2),allCTControlPoints_temp(:,3),'spline');
        allCTControlPoints = [allCTControlPoints; pt];
        allCTBifurcationPoints = [allCTBifurcationPoints; controlPoints(1).position'; controlPoints(size(controlPoints, 1)).position'];

    end





    % ptCTCloud_original = pointCloud(allCTControlPoints);
    % ptCTCloud = pctransform(ptCTCloud_original,tform_intial);

    ptCTCenterlineCloud = pointCloud(allCTControlPoints);
    ptCTBifurcationCloud = pointCloud(allCTBifurcationPoints);
    figure
    pcshow(ptCTCenterlineCloud);
    title('CT vessel feature');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold on
    plot3(allCTBifurcationPoints(:,1),allCTBifurcationPoints(:,2), allCTBifurcationPoints(:,3), '.r','MarkerSize',20)

    %%
    % Read vessel features (3DUS)
    % .json: all the vessel centerline sample points
    % .tsv: the propoerties of centerline. eg, vessel radius

    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\US_vessel_phantom_11\";
    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\Registration_results\20220613_1\";
    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\20220715_vessel_seg\VesselSegment3DUS_2\"
    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\20220715_vessel_seg\VesselSegment3DUS\"
    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\20220715_vessel_seg\Post_01\"
    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\20220715_vessel_seg\Post_02\"
    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\20220715_vessel_seg\Pre_01\"
    % dirPath = "E:\PROGRAM\Project_PhD\Slicing_verification\Vessel_based_registration\20220715_vessel_seg\Pre_03\"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HV02-02
    % 
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\registration\scene1\US_centerline\"
    % HV02-04
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US4\registration\scene1\US_centerline\';
    % HV02-03
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US3\registration\scene1\US_centerline\';
    % HV02-06
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US6\registration\Scene1\US_centerline\';
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\test\centerline\";
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\centerline\";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HV05-02 (RS)
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\registration\scene1\US_centerline\";
    % HV05-03
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\registration\scene1\US_centerline\';
    % HV05-06
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\registration\scene1\US_centerline\';
    % HV05-04
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_04\registration\scene1\US_centerline\';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% HV017-05
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_05\registration\scene1\US_centerline\";
% HV017-03
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_03\registration\scene1\US_centerline\";
% HV017-02
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_02\registration\scene1\US_centerline\";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% %PT013
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\US\registration\scene1\centerline\";

%PT012
dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\scene1\centerline\";

    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\initial\centerline\"
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\centerline_initial_1\"
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_04\centerline_initial\"
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\centerline_initial_1\"

    fileExtension = "*.json";
    searchingFiles = strcat(dirPath, fileExtension);
    jsonFiles = dir(searchingFiles);
    allUSControlPoints = [];
    allUSBifurcationPoints = [];
    % fileExtension = "*.tsv";
    % searchingFiles = strcat(dirPath, fileExtension);
    % tsvFiles = dir(searchingFiles);

    %
    % resample the vessel points for registration
    for i = 1:size(jsonFiles, 1)
        fileHeader = jsonFiles(i);
        fileName = fileHeader.name;
        fileFolder = fileHeader.folder;
        fileAbsoluteName = strcat(fileFolder, "\", fileName);
        fid = fopen(fileAbsoluteName); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        val = jsondecode(str);
        controlPoints = val.markups.controlPoints;
        allUSControlPoints_temp = [];
        for j = 1:size(controlPoints, 1)
            temp = controlPoints(j);
            allUSControlPoints_temp = [allUSControlPoints_temp; temp.position'];
    %         allUSControlPoints = [allUSControlPoints; temp.position'];
        end
    %     scaling = ceil(0.2*exp(allUSRadius(i)));
        scaling = 1;
        numOfInter = scaling*size(allUSControlPoints_temp,1);
        % interpolate using parametric splines
        pt = interparc(numOfInter,allUSControlPoints_temp(:,1),allUSControlPoints_temp(:,2),allUSControlPoints_temp(:,3),'spline');
        allUSControlPoints = [allUSControlPoints; pt];
        allUSBifurcationPoints = [allUSBifurcationPoints; controlPoints(1).position'; controlPoints(size(controlPoints, 1)).position'];
    end

    % flip points of 3D US
    allUSControlPoints = allUSControlPoints*[1, 0,  0; 0, 1, 0; 0, 0, 1];
    allUSBifurcationPoints = allUSBifurcationPoints* [1, 0,  0; 0, 1, 0; 0, 0, 1];
    ptUSCenterlineCloud = pointCloud(allUSControlPoints);
    ptUSBifurcationCloud = pointCloud(allUSBifurcationPoints);
    
    figure
    pcshow(ptUSCenterlineCloud);
    title('US vessel centerline feature');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    hold on
    plot3(allUSBifurcationPoints(:,1),allUSBifurcationPoints(:,2), allUSBifurcationPoints(:,3), '.b','MarkerSize',4)
    


    
 %%
 uniqueCTBifurcationPoints = unique(allCTBifurcationPoints(:,1));
 realCTBifurcationPoints = [];
for i = 1:size(uniqueCTBifurcationPoints)
    temp_index = find(allCTBifurcationPoints(:,1) == uniqueCTBifurcationPoints(i));
    temp_number = length(find(allCTBifurcationPoints(:,1) == uniqueCTBifurcationPoints(i)));
    if temp_number > 1
        realCTBifurcationPoints = [realCTBifurcationPoints; allCTBifurcationPoints(temp_index(1), :)];
    end
end

uniqueUSBifurcationPoints = unique(allUSBifurcationPoints(:,1));
 realUSBifurcationPoints = [];
for i = 1:size(uniqueUSBifurcationPoints)
    temp_index = find(allUSBifurcationPoints(:,1) == uniqueUSBifurcationPoints(i));
    temp_number = length(find(allUSBifurcationPoints(:,1) == uniqueUSBifurcationPoints(i)));
    if temp_number > 1
        realUSBifurcationPoints = [realUSBifurcationPoints; allUSBifurcationPoints(temp_index(1), :)];
    end
end

% % HV05-03
% pt1_US = [-62.4866, -22.7563, 155.5615];
% pt1_MR = [-82.8439, -34.7851, 148.8692];
% pt2_MR = [-46.2969, -39.5251, 172.9967];
% realCTBifurcationPoints = [realCTBifurcationPoints; [pt1_MR; pt2_MR]];
% realUSBifurcationPoints = [realUSBifurcationPoints; pt1_US];

% % % HV05-06
% pt1_MR = [-39.779, -3.207, 184.718];
% realCTBifurcationPoints = [realCTBifurcationPoints; pt1_MR];

% 
% HV0202
% pt1_MR = [-82.8736375, -34.200350000000018, 153.81295609756084; -74.31928576220443, 18.905829396717424, 142.42845350665164; 
%     -68.41848754882813, -1.059887409210205, 158.74697875976566;
%     -60.27059450001201, -10.300049372372677, 144.81264770984186;
%     -55.99279138875045, -12.298875182146972, 141.5804443359375;
%     -52.83482110974464, 18.71633147172726, 158.91799380106913;
%     -31.11661219788877, 12.74884966067053, 173.91796941082513]; % ->1
% realCTBifurcationPoints = [];
% realCTBifurcationPoints = [realCTBifurcationPoints; pt1_MR];


% % HV0204
% pt1_MR = [
%     -50.668796073045658, -15.951268145812307, 143.03692014219713;
%     -35.87393569946289, -44.61631011962891, 158.03689575195313;
%     -17.380360232484408, 6.703361801236298, 183.53685428853835;
%     -12.7569663657398, -8.553837959020925, 186.53684941048955];
% pt1_US = [
%     -48.96626068115235, -16.99935302734375, 142.9282054901123;
%     -33.37500808715821, -44.249608459472657, 160.01306354522704;
%     -20.837972948784939, 4.135195309777771, 185.81710029902517;
%     -13.311687636784939, -14.501320700888897, 185.81710029902517]
% realCTBifurcationPoints = [];
% realCTBifurcationPoints = [realCTBifurcationPoints; pt1_MR];
% realUSBifurcationPoints = [];
% realUSBifurcationPoints = [realUSBifurcationPoints; pt1_US];

% % % HV0203
% pt1_MR = [
%     -90.06019592285156, -40.51829586674491, 132.08247096803798;
%     -75.78275293273907, -24.12185464610002, 139.54427143697644;
%     -81.47172440488703, 15.315652714917519, 134.27211242437597;
%     -72.01378686115771, 19.475669942467748, 146.21630908739543;
%     -63.561432396496517, -34.7832916567622, 167.5143976312879;
%     -60.76371246188738, -11.163648536790234, 146.7025452947005];

% % realCTBifurcationPoints = [];
% realCTBifurcationPoints = [realCTBifurcationPoints; pt1_MR];

% HV02-06
% realCTBifurcationPoints = [];
% pt1_MR = [-67.54496851092847, -0.6271931539617981, 159.16308279234148;
%     -66.10380943129519, -7.690699330150721, 150.4608634563294;
%     -53.935984908974578, -28.069541767577918, 172.64047614494855;
%     -38.689104806536239, -20.436343533897224, 142.64052492543653;
%     -34.35033805343957, 13.885349457322463, 169.64048102299734]
% pt1_US = [-31.910297498960575, 11.963034688894809, 172.56719605436005];
% realUSBifurcationPoints = [realUSBifurcationPoints; pt1_US];
% realCTBifurcationPoints = [realCTBifurcationPoints; pt1_MR];

% HV017-05
% NONE


%     figure
%     pcshow(ptUSCenterlineCloud);
%     title('US vessel centerline feature');
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Z');
%     hold on
%     for i = 1: size(realUSBifurcationPoints, 1)
%         plot3(realUSBifurcationPoints(i,1),realUSBifurcationPoints(i,2), realUSBifurcationPoints(i,3), '.b','MarkerSize',4)
%         text(realUSBifurcationPoints(i,1),realUSBifurcationPoints(i,2), realUSBifurcationPoints(i,3), num2str(i), 'Color','w');
%     end
    
    
    
    
    figure (4)
    pcshowpair(ptCTCenterlineCloud,ptUSCenterlineCloud,'MarkerSize',20)
    hold on
    plot3(realCTBifurcationPoints(:,1),realCTBifurcationPoints(:,2), realCTBifurcationPoints(:,3), '.r','MarkerSize',20)
    plot3(realUSBifurcationPoints(:,1),realUSBifurcationPoints(:,2), realUSBifurcationPoints(:,3), '.b','MarkerSize',20)
    for i = 1: size(realUSBifurcationPoints, 1)
        plot3(realUSBifurcationPoints(i,1),realUSBifurcationPoints(i,2), realUSBifurcationPoints(i,3), '.b','MarkerSize',4)
        text(realUSBifurcationPoints(i,1),realUSBifurcationPoints(i,2), realUSBifurcationPoints(i,3), num2str(i), 'Color','w');
    end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds before registration')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')
 
 %%
 % HV0502
%  filter_US_bifuraction_index = [1];
%  % HV0503
%  filter_US_bifuraction_index = [5 6];
 
%   % HV0504
%  filter_US_bifuraction_index = [2 4];
   % HV0506
%  filter_US_bifuraction_index = [4];
    % HV0202
%  filter_US_bifuraction_index = [2, 3];
 
%   % HV0204
%    filter_US_bifuraction_index = [];
%      HV0203
%    filter_US_bifuraction_index = [3 10 13];
       % HV0206
%    filter_US_bifuraction_index = [6];

% % HV017-05
%    filter_US_bifuraction_index = [2 3 10 15];
% % HV017-03
%    filter_US_bifuraction_index = [2 3 4 5 15 16];

% HV017-02
   filter_US_bifuraction_index = [1 3 9 10 13 15 16 17 18 19];

% % PT013
% filter_US_bifuraction_index = [2]

% PT012
filter_US_bifuraction_index = [1 2 3 8 9]
   
 realUSBifurcationPoints_ = realUSBifurcationPoints;
 realUSBifurcationPoints_(filter_US_bifuraction_index,:) = [];
 matchedCTBifurcations = [];
 matchedBifurcationPointDistance = [];
 for i = 1: size(realUSBifurcationPoints_(:,1),1)
     delta = realUSBifurcationPoints_(i,:) - realCTBifurcationPoints;
     distance = sqrt(sum(delta.*delta, 2));
     [mini_distance, mini_index] = min(distance);
     matchedCTBifurcations = [matchedCTBifurcations; realCTBifurcationPoints(mini_index, :)];
     matchedBifurcationPointDistance = [matchedBifurcationPointDistance; mini_distance];
 end
 
xt = [matchedCTBifurcations(:,1), realUSBifurcationPoints_(:,1)];
yt = [matchedCTBifurcations(:,2), realUSBifurcationPoints_(:,2)];
zt = [matchedCTBifurcations(:,3), realUSBifurcationPoints_(:,3)];
    figure (3)
    pcshowpair(ptCTCenterlineCloud,ptUSCenterlineCloud,'MarkerSize',20)
    hold on
    plot3(matchedCTBifurcations(:,1),matchedCTBifurcations(:,2), matchedCTBifurcations(:,3), '.r','MarkerSize',20)
    plot3(realUSBifurcationPoints_(:,1),realUSBifurcationPoints_(:,2), realUSBifurcationPoints_(:,3), '.b','MarkerSize',20)
    for i = 1:size(realUSBifurcationPoints_,1)
     
     plot3(xt(i,:), yt(i,:), zt(i,:),'LineWidth', 1, 'Color', 'w');
     text(realUSBifurcationPoints_(i,1),realUSBifurcationPoints_(i,2), realUSBifurcationPoints_(i,3),num2str(matchedBifurcationPointDistance(i)), 'Color', 'w', 'FontSize', 20 )
    end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds before registration')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')
%%
figure
hist(matchedBifurcationPointDistance)
mean(matchedBifurcationPointDistance)
std(matchedBifurcationPointDistance)