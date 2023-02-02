close all; 
clear;
clc
addpath('..');

%% Step 1: Read vessels
studyDir = "E:\PROGRAM\Project_PhD\Registration\Results\2_vesselFeatureExtraction\LHV\"; % healthy volunteer study
caseName = "LHV-08";% healthy volunteer number
caseNumber = 1; % usually 4 3D ultrasound for registration
caseDir = strcat(studyDir, caseName);
listing = dir(strcat(caseDir,"\*.txt"));
casepath = strcat(listing(caseNumber).folder,'\', listing(caseNumber).name)

pathCell = readcell(casepath);
allCTSurfacePoints =  importdata(pathCell{3,2});
% 3DUS_01
allUSSurfacePoints = importdata(pathCell{6,2});
% 3DUS_02 
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_03\initial\vessel\vertices_vessel_US.txt');
        % 3DUS_03
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_02\initial\vessel\vertices_vessel_US.txt');
        % 3DUS_05
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_01\initial\vessel\vertices_vessel_US.txt');

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    initial_transform_slicer = load(pathCell{8,2});
    initial_transform_invert = [reshape(initial_transform_slicer.AffineTransform_double_3_3(1:9),[3,3])' initial_transform_slicer.AffineTransform_double_3_3(10:12); 0 0 0 1];
    initial_transform = inv(initial_transform_invert);
    rotation_initial = initial_transform(1:3, 1:3);
    translation_initial = initial_transform(1:3, 4);
    
    allUSSurfacePoints_inital = (rotation_initial*allUSSurfacePoints' + translation_initial)';
    allCTSurfacePoints_inital = allCTSurfacePoints;
    
  
    figure(1)
    plot3(allCTSurfacePoints_inital(:,1), allCTSurfacePoints_inital(:,2), allCTSurfacePoints_inital(:,3),'.b','MarkerSize',3); 
    
    hold on;
%     figure
    plot3(allUSSurfacePoints_inital(:,1), allUSSurfacePoints_inital(:,2), allUSSurfacePoints_inital(:,3),'.r','MarkerSize',3); 
    daspect([1 1 1]); 
    grid on;
    title('Original vessel surface point sets in CT/MRI and US','FontSize',18);
    lgd = legend('MRI','US');
    lgd.FontSize = 12;
%     lgd.FontWeight = 'bold'
    

%%
    % filter centerlines

    % compute the FOV of moving images (3D US)
    x_up = max(allUSSurfacePoints_inital(:,1));
    x_down = min(allUSSurfacePoints_inital(:,1));
    y_up = max(allUSSurfacePoints_inital(:,2));
    y_down = min(allUSSurfacePoints_inital(:,2));
    z_up = max(allUSSurfacePoints_inital(:,3));
    z_down = min(allUSSurfacePoints_inital(:,3));
    center = [(x_up - x_down)/2 + x_down, (y_up - y_down)/2 + y_down, (z_up - z_down)/2 + z_down];

    % Customize the FOV for fixed images (CT, MRI)
    scale = 1.2;
    x_up_new = center(1) + (x_up - x_down)*scale/2;
    x_down_new = center(1) - (x_up - x_down)*scale/2;
    y_up_new = center(2) + (y_up - y_down)*scale/2;
    y_down_new = center(2) - (y_up - y_down)*scale/2;
    z_up_new = center(3) + (z_up - z_down)*scale/2;
    z_down_new = center(3) - (z_up - z_down)*scale/2;

    % filter the CT/MRI
    allCTSurfacePoints_FOV = allCTSurfacePoints_inital;
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) > x_up_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) < x_down_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];

%     ptUSSurfaceCloud = pointCloud(allCTSurfacePoints);
    ptCTSurfaceCloud_FOV = pointCloud(allCTSurfacePoints_FOV);
    
    figure(2)
    plot3(allCTSurfacePoints_FOV(:,1), allCTSurfacePoints_FOV(:,2), allCTSurfacePoints_FOV(:,3),'.b','MarkerSize',4); 
    hold on;
    plot3(allUSSurfacePoints_inital(:,1), allUSSurfacePoints_inital(:,2), allUSSurfacePoints_inital(:,3),'.r','MarkerSize',4); 
    daspect([1 1 1]); 
    grid on;
    title('Point clouds before registration after FOV','FontSize',18);
    lgd = legend('MRI','US');
    lgd.FontSize = 12;


    x   =sprintf(pathCell{9,2},pwd);
    y   =sprintf(pathCell{10,2},pwd);

    writematrix(allCTSurfacePoints_FOV,x ,'Delimiter','tab')
    writematrix(allUSSurfacePoints_inital,y,'Delimiter','tab')
%%
%     % input files
studyDir = "E:\PROGRAM\Project_PhD\Registration\Results\2_vesselFeatureExtraction\LHV\"; % healthy volunteer study
caseName = "LHV-08";% healthy volunteer number
caseNumber = 1; % usually 4 3D ultrasound for registration
caseDir = strcat(studyDir, caseName);
listing = dir(strcat(caseDir,"\*.txt"));
casepath = strcat(listing(caseNumber).folder,'\', listing(caseNumber).name)

pathCell = readcell(casepath);
    x   =sprintf(pathCell{9,2},pwd);
    y   =sprintf(pathCell{10,2},pwd);


    fnm =sprintf('%s/bcpd_package/bcpd',                pwd);
    fnw =sprintf('%s/bcpd_package/win/bcpd.exe',        pwd);
    if(ispc) bcpd=fnw; else bcpd=fnm; end;

    %
    % parameters
    omg ='0.5'; % outlier
    bet ='2.0';
    lmd ='1e9';
    gma ='0.1';

    K   ='70'; % 
    J   ='300';
    f   ='0.3';
    c   ='1e-9';% itration tolerance
    n   ='1000'; % maximum iteration number
    N = '300'; % minum iteration number
    nrm ='y'; % normalization
    dwnx ='x,15000,0.02'; % downsampling 
    dwny ='y,15000,0.02'; % downsampling
    % 
    outputName = 'Output';
    
    % save out settings
cpd_setting_cell = {
        'omg', omg;
        'bet', bet;
        'lmd', lmd;
        'gma', gma;
        'K', K;
        'J', J;
        'f', f;
        'c', c;
        'n', n;
        'N', N;
        'nrm', nrm;
        'dwnx', dwnx;
        'dwny', dwny;
        
        }

    writecell(cpd_setting_cell,pathCell{13,2})

    % execution
    prm1=sprintf('-w%s -b%s -l%s -g%s',omg,bet,lmd,gma);
%     prm2=sprintf('-J%s -K%s -p -f%s -u%s -D%s -D%s',J,K,f,nrm,dwnx,dwny);
    prm2=sprintf('-J%s -K%s -p -d5 -e0.3 -f%s -u%s',J,K,f,nrm);
    prm3=sprintf('-c%s -N%s -n%s -h -r1',c,N,n);
    cmd =sprintf('%s -x%s -y%s %s %s %s -sYT',bcpd,x,y,prm1,prm2,prm3);
    system(cmd); 
    optpath3;



    
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
studyDir = "E:\PROGRAM\Project_PhD\Registration\Results\2_vesselFeatureExtraction\LHV\"; % healthy volunteer study
caseName = "LHV-08";% healthy volunteer number
caseNumber = 4; % usually 4 3D ultrasound for registration
caseDir = strcat(studyDir, caseName);
listing = dir(strcat(caseDir,"\*.txt"));
casepath = strcat(listing(caseNumber).folder,'\', listing(caseNumber).name)

pathCell = readcell(casepath); 
allCTSurfacePoints =  importdata(pathCell{3,2});
% 3DUS_01
allUSSurfacePoints = importdata(pathCell{6,2});
initial_transform_slicer = load(pathCell{8,2});
initial_transform_invert = [reshape(initial_transform_slicer.AffineTransform_double_3_3(1:9),[3,3])' initial_transform_slicer.AffineTransform_double_3_3(10:12); 0 0 0 1];
initial_transform = inv(initial_transform_invert);
rotation_initial = initial_transform(1:3, 1:3);
translation_initial = initial_transform(1:3, 4);

allUSSurfacePoints_inital = (rotation_initial*allUSSurfacePoints' + translation_initial)';
allCTSurfacePoints_inital = allCTSurfacePoints;
    
%     normX = sprintf('%s/output_normX.txt',pwd);
%     normY = sprintf('%s/output_normY.txt',pwd);
%     output_y = sprintf('%s/output_y.txt',pwd);
%     load(normX)
%     load(normY)
%     load(output_y)
   % compute the FOV of moving images (3D US)
    x_up = max(allUSSurfacePoints_inital(:,1));
    x_down = min(allUSSurfacePoints_inital(:,1));
    y_up = max(allUSSurfacePoints_inital(:,2));
    y_down = min(allUSSurfacePoints_inital(:,2));
    z_up = max(allUSSurfacePoints_inital(:,3));
    z_down = min(allUSSurfacePoints_inital(:,3));
    center = [(x_up - x_down)/2 + x_down, (y_up - y_down)/2 + y_down, (z_up - z_down)/2 + z_down];

    % Customize the FOV for fixed images (CT, MRI)
    scale = 1.2;
    
    
    x_up_new = center(1) + (x_up - x_down)*scale/2;
    x_down_new = center(1) - (x_up - x_down)*scale/2;
    y_up_new = center(2) + (y_up - y_down)*scale/2;
    y_down_new = center(2) - (y_up - y_down)*scale/2;
    z_up_new = center(3) + (z_up - z_down)*scale/2;
    z_down_new = center(3) - (z_up - z_down)*scale/2;

    % filter the CT/MRI
    allCTSurfacePoints_FOV = allCTSurfacePoints_inital;
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) > x_up_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) < x_down_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];
   

%     [normMR, normMatrixMR, muX, sgmX] = normalization(allCTSurfacePoints, 1, allUSSurfacePoints);
%     [normUS, normMatrixUS, muY, sgmY] = normalization(allUSSurfacePoints, 1, allUSSurfacePoints);

    [normMR, normMatrixMR, muX, sgmX] = normalization(allCTSurfacePoints_FOV, 1, allUSSurfacePoints_inital);
    [normUS, normMatrixUS, muY, sgmY] = normalization(allUSSurfacePoints_inital, 1, allUSSurfacePoints_inital);
%     figure(4)
%     scatter3(output_y(:,1), output_y(:,2), output_y(:,3),'r');
%     hold on
%     scatter3(MR_pt(:,1), MR_pt(:,2), MR_pt(:,3),'g');
%     hold off



    %
    R = sprintf('%s/output_R.txt',pwd);
    t = sprintf('%s/output_t.txt',pwd);
    s = sprintf('%s/output_s.txt',pwd);
    load(R); load(t); load(s);

    % apply transformation

    s0 = output_s*(sgmX/sgmY)
    R0 = output_R
%     t0 = sgmX*output_t + muX' - output_s*(sgmX/sgmY)*output_R*muY'
    t0 = sgmX*output_t + muX' - (sgmX/sgmY)*output_R*muY'
    
%     movingImg = (R0*s0*allUSSurfacePoints_inital' + t0)';
    movingImg = (R0*allUSSurfacePoints_inital' + t0)';


%     figure(5)
%     scatter3(movingImg(:,1), movingImg(:,2), movingImg(:,3),'r');
%     hold on
%     scatter3(MR_pt(:,1), MR_pt(:,2), MR_pt(:,3),'g');
%     hold off

%     f3=figure('Name','Before/After Registration','NumberTitle','off'); set(f2,'Position',w1);
%     subplot(1,2,1);
%     plot3(output_y(:,1), output_y(:,2), output_y(:,3),'.r','MarkerSize',3); hold on;
%     plot3(MR_pt(:,1), MR_pt(:,2), MR_pt(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
%     title('Before Registration','FontSize',18);
%     subplot(1,2,2);
%     plot3(output_y(:,1), output_y(:,2), output_y(:,3),'.r','MarkerSize',3); hold on;
%     plot3(movingImg(:,1), movingImg(:,2), movingImg(:,3),'.b','MarkerSize',3); daspect([1 1 1]); grid on;
%     title('After Registration','FontSize',18);


     % compute the FOV of moving images (3D US)
    x_up = max(allUSSurfacePoints_inital(:,1));
    x_down = min(allUSSurfacePoints_inital(:,1));
    y_up = max(allUSSurfacePoints_inital(:,2));
    y_down = min(allUSSurfacePoints_inital(:,2));
    z_up = max(allUSSurfacePoints_inital(:,3));
    z_down = min(allUSSurfacePoints_inital(:,3));
    center = [(x_up - x_down)/2 + x_down, (y_up - y_down)/2 + y_down, (z_up - z_down)/2 + z_down];

    % Customize the FOV for fixed images (CT, MRI)
    scale = 2;
    x_up_new = center(1) + (x_up - x_down)*scale/2;
    x_down_new = center(1) - (x_up - x_down)*scale/2;
    y_up_new = center(2) + (y_up - y_down)*scale/2;
    y_down_new = center(2) - (y_up - y_down)*scale/2;
    z_up_new = center(3) + (z_up - z_down)*scale/2;
    z_down_new = center(3) - (z_up - z_down)*scale/2;

    % filter the CT/MRI
    allCTSurfacePoints_FOV = allCTSurfacePoints_inital;
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) > x_up_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) < x_down_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];
    
    f3=figure('Name','Before/After Coarse Registration','NumberTitle','off'); set(f3,'Position',w1);
    subplot(1,2,1);
    plot3(allUSSurfacePoints_inital(:,1),allUSSurfacePoints_inital(:,2),allUSSurfacePoints_inital(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTSurfacePoints_FOV(:,1),allCTSurfacePoints_FOV(:,2),allCTSurfacePoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('Before Coarse Registration','FontSize',18);
    subplot(1,2,2);
    plot3(movingImg(:,1),movingImg(:,2),movingImg(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTSurfacePoints_FOV(:,1),allCTSurfacePoints_FOV(:,2),allCTSurfacePoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('After Coarse Registration','FontSize',18);
    
    vessel_transformationCPD_LPS = [R0 t0; 0 0 0 1];
    
%% Step 2: Read vessel centerlines

    dirPathMR = pathCell{4,2};
    fileExtension = "*.json";
    searchingFiles = strcat(dirPathMR, fileExtension);
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
    scatter3(allCTBifurcationPoints(:,1),allCTBifurcationPoints(:,2), allCTBifurcationPoints(:,3), 'LineWidth', 1)

    %%
    % Read vessel features (3DUS)
    
    dirPathUS = pathCell{7,2};
    fileExtension = "*.json";
    searchingFiles = strcat(dirPathUS, fileExtension);
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
    ptUSCloud = pointCloud(allUSControlPoints);
    ptUSBifurcationCloud = pointCloud(allUSBifurcationPoints);
    
    % initial transform US control points from CPD surface registration
    allUSControlPoints_inital = (rotation_initial*allUSControlPoints' + translation_initial)'; 
    %ptUSCenterlineCloud_initial = pointCloud(allUSControlPoints_initial);
    
    
    allUSControlPoints_ = (vessel_transformationCPD_LPS(1:3, 1:3)*allUSControlPoints_inital' + vessel_transformationCPD_LPS(1:3,4))'; 
    ptUSCenterlineCloud_ = pointCloud(allUSControlPoints_);
    
    
    figure 
    pcshowpair(ptCTCenterlineCloud,ptUSCenterlineCloud_,'MarkerSize',50)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds before registration')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')
%%
    % filter centerlines

    % compute the FOV of moving images (3D US)
    x_up = max(allUSControlPoints_(:,1));
    x_down = min(allUSControlPoints_(:,1));
    y_up = max(allUSControlPoints_(:,2));
    y_down = min(allUSControlPoints_(:,2));
    z_up = max(allUSControlPoints_(:,3));
    z_down = min(allUSControlPoints_(:,3));
    center = [(x_up - x_down)/2 + x_down, (y_up - y_down)/2 + y_down, (z_up - z_down)/2 + z_down];

    % Customize the FOV for fixed images (CT, MRI)
    scale = 1.1;
    x_up_new = center(1) + (x_up - x_down)*scale/2;
    x_down_new = center(1) - (x_up - x_down)*scale/2;
    y_up_new = center(2) + (y_up - y_down)*scale/2;
    y_down_new = center(2) - (y_up - y_down)*scale/2;
    z_up_new = center(3) + (z_up - z_down)*scale/2;
    z_down_new = center(3) - (z_up - z_down)*scale/2;

    % filter the CT/MRI
    allCTControlPoints_FOV = allCTControlPoints;
    allCTControlPoints_FOV(find(allCTControlPoints_FOV(:,1) > x_up_new), :) = [];
    allCTControlPoints_FOV(find(allCTControlPoints_FOV(:,1) < x_down_new), :) = [];
    allCTControlPoints_FOV(find(allCTControlPoints_FOV(:,2) > y_up_new),:) = [];
    allCTControlPoints_FOV(find(allCTControlPoints_FOV(:,2) < y_down_new),:) = [];
    allCTControlPoints_FOV(find(allCTControlPoints_FOV(:,2) > y_up_new),:) = [];
    allCTControlPoints_FOV(find(allCTControlPoints_FOV(:,2) < y_down_new),:) = [];


%     ptUSCloud = pointCloud(allUSControlPoints);
    ptCTCenterlineCloud_FOV = pointCloud(allCTControlPoints_FOV);

    figure (5)
    pcshowpair(ptCTCenterlineCloud_FOV,ptUSCenterlineCloud_,'MarkerSize',50)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds before registration after FOV')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')
    
    ptCTCenterlineCloudSampled = pcdownsample(ptCTCenterlineCloud_FOV, 'random', 1);
    % ptCTCloudSampled = pcdownsample(ptCTCloud, 'random', 1);

    ptUSCenterlineCloudSampled = pcdownsample(ptUSCenterlineCloud_, 'random', 1);

    figure(6)
    pcshowpair(ptCTCenterlineCloud_FOV,ptUSCenterlineCloud_,'MarkerSize',50)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds before registration')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')

    figure(5)
    pcshowpair(ptCTCenterlineCloudSampled,ptUSCenterlineCloudSampled,'MarkerSize',50)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds after sampling and FOV' )
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')

 %% ICP

    N_iter = 1000;
    InlierRatio = 0.5; 
    Tolerance = [0.001, 0.005];
   
    tic;
    [tform, movingReg, rmse] = pcregistericp(ptUSCenterlineCloudSampled,ptCTCenterlineCloudSampled, 'MaxIterations', N_iter, 'InlierRatio', InlierRatio, 'Tolerance', Tolerance, 'Verbos', true);
    toc;
    movingReg = pctransform(ptUSCenterlineCloudSampled,tform);
    % ptUSCloud = pointCloud(allUSControlPoints);
    figure (6)
    pcshowpair(ptCTCenterlineCloudSampled,movingReg,'MarkerSize',50)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds after registration (ICP)')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')

    disp('RMSE:')
    rmse
    disp('Rigid transform(LPS): ')
    centerline_transformationICP_LPS = tform.T'
    disp('Rigid transform(RAS): ')
    transform_RAS = [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1]
    centerline_registrationTransform_RAS = transform_RAS*tform.T'*transform_RAS
%     % test
%     movingTest = centerline_transformationICP_LPS(1:3, 1:3)*allUSControlPoints_' + centerline_transformationICP_LPS(1:3,4); 
%     figure(8)
%     ptcloudTest = pointCloud(movingTest');
%     pcshowpair(ptcloudTest,movingReg,'MarkerSize',50)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     title('Point clouds after registration')
%     legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
%     legend('Location','southoutside')

%%

    % save out settings
icp_setting_cell = {
        'MaxIterations', num2str(N_iter);
        'InlierRatio', num2str(InlierRatio);
        'Tolerance', strcat('[', num2str(Tolerance(1)), ',', num2str(Tolerance(2)), ']');
        }

    writecell(icp_setting_cell,pathCell{14,2})
%% CPD
    N_iter = 1000;
    tic;
    [tform, movingReg, rmse] = pcregistercpd(ptUSCenterlineCloudSampled,ptCTCenterlineCloudSampled, 'Transform', 'Rigid', 'MaxIterations', N_iter, 'OutlierRatio', 0.2, 'Tolerance', 1e-8, 'Verbos', true);
    toc;
    movingReg = pctransform(ptUSCenterlineCloudSampled,tform);
    % ptUSCloud = pointCloud(allUSControlPoints);
    figure (7)
    pcshowpair(movingReg,ptCTCenterlineCloudSampled,'MarkerSize',50)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds after registration')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')

    disp('RMSE:')
    rmse
    disp('Rigid transform(LPS): ')
    centerline_transformationCPD_LPS = tform.T'

    disp('Rigid transform(RAS): ')
    transform_RAS = [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1]
    centerline_registrationTransform_RAS = transform_RAS*tform.T'*transform_RAS



%% compute the final transformation
    R = centerline_transformationICP_LPS(1:3, 1:3)*vessel_transformationCPD_LPS(1:3, 1:3);
    t = centerline_transformationICP_LPS(1:3, 1:3)*vessel_transformationCPD_LPS(1:3, 4) + centerline_transformationICP_LPS(1:3, 4);
    registrationMatrix = [R, t; 0 0 0 1]
%     registrationMatrix = (centerline_transformationICP_LPS'*vessel_transformationCPD_LPS')';
    movingTest = (registrationMatrix(1:3, 1:3)*allUSControlPoints_inital' + registrationMatrix(1:3,4))'; 
%     movingTest = R*allUSControlPoints' + t; 
    
%     figure(8)
%     ptcloudTest = pointCloud(movingTest);
%     pcshowpair(ptcloudTest,movingReg,'MarkerSize',50)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     title('Point clouds after registration')
%     legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
%     legend('Location','southoutside')
%     
    f2=figure('Name','Before/After Fine Registration','NumberTitle','off'); set(f2,'Position',w1);
    subplot(1,2,1);
    plot3(allUSControlPoints_(:,1),allUSControlPoints_(:,2),allUSControlPoints_(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTControlPoints_FOV(:,1),allCTControlPoints_FOV(:,2),allCTControlPoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('Before Fine Registration','FontSize',18);
    subplot(1,2,2);
    plot3(movingTest(:,1),movingTest(:,2),movingTest(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTControlPoints_FOV(:,1),allCTControlPoints_FOV(:,2),allCTControlPoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('After Fine Registration','FontSize',18);

% save out the file and settings
%     transform_RAS = [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1]
%     registrationMatrix_RAS = transform_RAS*registrationMatrix*transform_RAS;
%     
%     save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\combined_transformation_ini_LPS.txt', 'registrationMatrix', '-ascii');
%     save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\combined_transformation_ini_RAS.txt', 'registrationMatrix_RAS', '-ascii');
%     
%     save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\vessel_transformation_ini_LPS.txt', 'vessel_transformationCPD_LPS', '-ascii');
%     save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\centerline_transformationCPD_ini_LPS.txt', 'centerline_transformationCPD_LPS', '-ascii');

    SlicerTransformationConverter(vessel_transformationCPD_LPS, pathCell{11,2})
    SlicerTransformationConverter(registrationMatrix, pathCell{12,2})
    
    %     save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\registration\transformation_ini_RAS.txt', 'registrationMatrix_RAS', '-ascii');





