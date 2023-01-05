close all; 
clear;
clc
addpath('..');
% addpath('E:\PROGRAM\Project_PhD\Registration\bcpd\bcpd-master\demo')

%% Step 1: Read vessels

    % vessel surface registration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HV05
%
    %    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\Vessel\vertices_vessel_MR.txt');
    %    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\Vessel\vertices_vessel_IVC_MR.txt');
    %    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\Vessel\vertices_vessel_MR_refined.txt');

    %     % HV05-02 (RS)
    %    US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\initial\Vessel\vertices_vessel_US.txt');
    %       US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\initial\Vessel\forPaper\vertices_vessel_US.txt');

        % HV05-03
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\initital\Vessel\vertices_vessel_US.txt');
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\initital\Vessel\vertices_vessel_US_trimmed.txt');
        % HV05-06
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\initial_correct\Vessel_1\vertices_vessel_US.txt');

    %     HV05-04
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_04\initial\Vessel\vertices_vessel_US.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%
%HV02
%
    %     MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Vessel\vertices_vessel_MR.txt');
    %     MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Vessel\vertices_vessel_MR_IVC.txt');
    %       MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Vessel\vertices_vessel_MR_refined_1.txt');
    %       MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Vessel\vertices_vessel_MR_IVC_refined_1.txt');
    %     % HV02-02
    %       US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\initial\Vessel\vertices_vessel_US_refined.txt');
    %       
          % HV02-04
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US4\initial\Vessel\vertices_vessel_US.txt');
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US4\initial\Vessel_manual\vertices_vessel_US.txt');  
          % HV02-03
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US3\initial\Vessel_manual\vertices_vessel_US.txt');
        % HV02-06
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US6\initial\Vessel\vertices_vessel_US1.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HV009
%
pathCell = readcell("E:\PROGRAM\Project_PhD\Registration\Results\2_vesselFeatureExtraction\LHV\LHV-09\3DUS_1\LHV09_3DUS1.txt");

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
%
% HV017
%
    %    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR\Vessel\vertices_vessel_MR.txt');

       % HV017-05
    %    US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_05\initial\vessel\vertices_vessel_US_refined.txt');
       % HV017-03   
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_03\initial\vessel\vertices_vessel_US.txt');
        % HV017-02
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_02\initial\vessel\vertices_vessel_US.txt');
        % HV017-01
    %     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_01\initial\vessel\vertices_vessel_US.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%PT013
%
    %    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\MR\vessel\vertices_vessel_MR.txt');
    %    US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\US\initial\vessel\vertices_vessel_US.txt');
%
% PT012
%
%    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\CT\vessel\vertices_vessel_CT.txt');
%    US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\vessel\vertices_vessel_US.txt');
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
%%    

%     registrationMatrix = [0.381902 0.909698 -0.163095 19.2727; 
%     -0.286736 -0.0511372 -0.956644 -15.1374; 
%     -0.878598 0.412109 0.241313 8.24036; 
%     0 0 0 1];
%     transform_RAS = [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1];
%     registrationMatrix_RAS = transform_RAS*registrationMatrix*transform_RAS
%     inv(registrationMatrix_RAS)

    initial_transform_slicer = load(pathCell{8,2});
    initial_transform_invert = [reshape(initial_transform_slicer.AffineTransform_double_3_3(1:9),[3,3])' initial_transform_slicer.AffineTransform_double_3_3(10:12); 0 0 0 1];
    initial_transform = inv(initial_transform_invert);
    rotation_initial = initial_transform(1:3, 1:3);
    translation_initial = initial_transform(1:3, 4);
    
    allUSSurfacePoints_inital = (rotation_initial*allUSSurfacePoints' + translation_initial)';
    allCTSurfacePoints_inital = allCTSurfacePoints;
    
  
    ptCTSurfaceCloud = pointCloud(allCTSurfacePoints_inital);
    ptUSSurfaceCloud = pointCloud(allUSSurfacePoints_inital);
    
%     figure (1)
%     pcshowpair(ptCTSurfaceCloud,ptUSSurfaceCloud,'MarkerSize',50)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     title('Point clouds before registration before FOV')
%     legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')

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
    scale = 1.1;
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
    plot3(allCTSurfacePoints_FOV(:,1), allCTSurfacePoints_FOV(:,2), allCTSurfacePoints_FOV(:,3),'.r','MarkerSize',4); 
    hold on;
    plot3(allUSSurfacePoints_inital(:,1), allUSSurfacePoints_inital(:,2), allUSSurfacePoints_inital(:,3),'.b','MarkerSize',4); 
    daspect([1 1 1]); 
    grid on;
    title('Point clouds before registration after FOV','FontSize',18);
    lgd = legend('MRI','US');
    lgd.FontSize = 12;
    
%     figure (2)
%     pcshowpair(ptCTSurfaceCloud_FOV,ptUSSurfaceCloud,'MarkerSize',50)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     title('Point clouds before registration after FOV')
%     legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
%     legend('Location','southoutside')

    x   =sprintf(pathCell{9,2},pwd);
    y   =sprintf(pathCell{10,2},pwd);

    writematrix(allCTSurfacePoints_FOV,x ,'Delimiter','tab')
    writematrix(allUSSurfacePoints_inital,y,'Delimiter','tab')
%%
%     % input files
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
    % HV05
%     MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\Vessel\vertices_vessel_MR.txt');
%     MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\Vessel\vertices_vessel_IVC_MR.txt');
%     MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\Vessel\vertices_vessel_MR_refined.txt');
%     % HV05-02
%     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\initial\Vessel\vertices_vessel_US.txt');
%     % HV05-03
%     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\initital\Vessel\vertices_vessel_US_trimmed.txt');
%     % HV05-06
%     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\initial_correct\Vessel_1\vertices_vessel_US.txt');
     % HV05-04
%     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_04\initial\Vessel\vertices_vessel_US.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%     % HV02
%     MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Vessel\vertices_vessel_MR.txt');
%       MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Vessel\vertices_vessel_MR_IVC.txt');
%       MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Vessel\vertices_vessel_MR_refined_1.txt');
%     % HV02-02
%       US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\initial\Vessel\vertices_vessel_US_refined.txt');
      % HV02-04
%     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US4\initial\Vessel_manual\vertices_vessel_US.txt'); 
      % HV02-03
%       US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US3\initial\Vessel_manual\vertices_vessel_US.txt');

  % HV02-06
%     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US6\initial\Vessel\vertices_vessel_US1.txt');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% % PT013
%    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\MR\vessel\vertices_vessel_MR.txt');
%    US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\US\initial\vessel\vertices_vessel_US.txt');
%  
% % PT012
%   MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\CT\vessel\vertices_vessel_CT.txt');
%   US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\vessel\vertices_vessel_US.txt');

   
    %
    % HV017
    %
%    MR_vessel =  importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR\Vessel\vertices_vessel_MR.txt');
   
   % HV017-05
%    US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_05\initial\vessel\vertices_vessel_US_refined.txt');
    % HV017-03
%    US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_03\initial\vessel\vertices_vessel_US.txt');
    % HV017-02
%     US_vessel = importdata('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_02\initial\vessel\vertices_vessel_US.txt');
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
 pathCell = readcell("E:\PROGRAM\Project_PhD\Registration\Results\2_vesselFeatureExtraction\LHV\LHV-09\3DUS_1\LHV09_3DUS1.txt");

       allCTSurfacePoints =  importdata(pathCell{3,2});
       % 3DUS_01
       allUSSurfacePoints = importdata(pathCell{6,2});
 
 
    %     normX = sprintf('%s/output_normX.txt',pwd);
%     normY = sprintf('%s/output_normY.txt',pwd);
%     output_y = sprintf('%s/output_y.txt',pwd);
%     load(normX)
%     load(normY)
%     load(output_y)
    

    [normMR, normMatrixMR, muX, sgmX] = normalization(allCTSurfacePoints, 1, allUSSurfacePoints);
    [normUS, normMatrixUS, muY, sgmY] = normalization(allUSSurfacePoints, 1, allUSSurfacePoints);

%     figure(4)
%     scatter3(output_y(:,1), output_y(:,2), output_y(:,3),'r');
%     hold on
%     scatter3(MR_pt(:,1), MR_pt(:,2), MR_pt(:,3),'g');
%     hold off


    % back up

%     ptUSCenterlineCloudSampled = pointCloud(US_2);
%     ptCTCenterlineCloudSampled = pointCloud(output_y);
%     N_iter = 1000;
%     [tform, movingReg, rmse] = pcregistercpd(ptUSCenterlineCloudSampled,ptCTCenterlineCloudSampled, 'Transform', 'Rigid', 'MaxIterations', N_iter, 'OutlierRatio', 0.02, 'Tolerance', 1e-8, 'Verbos', true);
% 
%     movingReg = pctransform(ptUSCenterlineCloudSampled,tform);
%     % ptUSCloud = pointCloud(allUSControlPoints);
%     figure (6)
%     pcshowpair(movingReg,ptCTCenterlineCloudSampled,'MarkerSize',50)
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     title('Point clouds after registration')
%     legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
%     legend('Location','southoutside')
% 
%     disp('RMSE:')
%     rmse
%     disp('Rigid transform(LPS): ')
%     centerline_transformationCPD_LPS = tform.T'
%     % figure
%     % scatter3(output_y(:,1), output_y(:,2), output_y(:,3),'r');
%     % 
%     % figure
%     % scatter3(US_2(:,1), US_2(:,2), US_2(:,3),'r');
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

%     movingImg = (R0*s0*US_vessel' + t0)';
    movingImg = (R0*US_vessel' + t0)';


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
    x_up = max(allUSSurfacePoints(:,1));
    x_down = min(allUSSurfacePoints(:,1));
    y_up = max(allUSSurfacePoints(:,2));
    y_down = min(allUSSurfacePoints(:,2));
    z_up = max(allUSSurfacePoints(:,3));
    z_down = min(allUSSurfacePoints(:,3));
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
    allCTSurfacePoints_FOV = allCTSurfacePoints;
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) > x_up_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,1) < x_down_new), :) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) > y_up_new),:) = [];
    allCTSurfacePoints_FOV(find(allCTSurfacePoints_FOV(:,2) < y_down_new),:) = [];
    
    f3=figure('Name','Before/After Coarse Registration','NumberTitle','off'); set(f3,'Position',w1);
    subplot(1,2,1);
    plot3(allUSSurfacePoints(:,1),allUSSurfacePoints(:,2),allUSSurfacePoints(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTSurfacePoints_FOV(:,1),allCTSurfacePoints_FOV(:,2),allCTSurfacePoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('Before Coarse Registration','FontSize',18);
    subplot(1,2,2);
    plot3(movingImg(:,1),movingImg(:,2),movingImg(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTSurfacePoints_FOV(:,1),allCTSurfacePoints_FOV(:,2),allCTSurfacePoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('After Coarse Registration','FontSize',18);
    
    vessel_transformationCPD_LPS = [R0 t0; 0 0 0 1];
    
%% Step 2: Read vessel centerlines
    % HV-05
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\centerline\"
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR\centerline_refined\"
    
    % HV-02
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR\Centerline\";
% % HV17
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR\Centerline\";
% PT13
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\MR\centerline\";

% % PT12
% dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\CT\centerline\"

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
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\initial\Centerline\"
    % HV02-04
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US4\initial\Centerline_manual\';
    % HV02-03
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US3\initial\Centerline_manual\';
    % HV02-06
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US6\initial\Centerline\';
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\test\centerline\";
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-02\MR_US2\centerline\";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % HV05-02 (RS)
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\initial\centerline\";
    % HV05-03
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\initital\Centerline\';
    % HV05-06
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\initial_correct\Centerline_1\';
    % HV05-04
%     dirPath = 'E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_04\initial\Centerline\';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % HV17-05 (RS)
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_05\initial\centerline\";
    % HV17-03
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_03\initial\centerline\";
    % HV17-02
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-17\MR_US_02\initial\centerline\";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% PT13
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-13\US\initial\centerline\";
% PT12
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\centerline\"
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\initial\centerline\"
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_03\centerline_initial_1\"
    % dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_04\centerline_initial\"
%     dirPath = "E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_06\centerline_initial_1\"
    
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
    
    allUSControlPoints_ = (vessel_transformationCPD_LPS(1:3, 1:3)*allUSControlPoints' + vessel_transformationCPD_LPS(1:3,4))'; 
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
    scale = 1;
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
    %%
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
    tic;
    [tform, movingReg, rmse] = pcregistericp(ptUSCenterlineCloudSampled,ptCTCenterlineCloudSampled, 'MaxIterations', N_iter, 'InlierRatio', 0.5, 'Tolerance', [0.001, 0.005], 'Verbos', true);
    toc;
    movingReg = pctransform(ptUSCenterlineCloudSampled,tform);
    % ptUSCloud = pointCloud(allUSControlPoints);
    figure (6)
    pcshowpair(movingReg,ptCTCenterlineCloudSampled,'MarkerSize',50)
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

%% CPD
    N_iter = 1000;
    tic;
    [tform, movingReg, rmse] = pcregistercpd(ptUSCenterlineCloudSampled,ptCTCenterlineCloudSampled, 'Transform', 'Rigid', 'MaxIterations', N_iter, 'OutlierRatio', 0, 'Tolerance', 1e-8, 'Verbos', true);
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
    movingTest = (registrationMatrix(1:3, 1:3)*allUSControlPoints' + registrationMatrix(1:3,4))'; 
%     movingTest = R*allUSControlPoints' + t; 
    
    figure(8)
    ptcloudTest = pointCloud(movingTest);
    pcshowpair(ptcloudTest,movingReg,'MarkerSize',50)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Point clouds after registration')
    legend({'CT centerline point cloud','US centerline point cloud'},'TextColor','w')
    legend('Location','southoutside')
    
    f2=figure('Name','Before/After Fine Registration','NumberTitle','off'); set(f2,'Position',w1);
    subplot(1,2,1);
    plot3(allUSControlPoints_(:,1),allUSControlPoints_(:,2),allUSControlPoints_(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTControlPoints_FOV(:,1),allCTControlPoints_FOV(:,2),allCTControlPoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('Before Fine Registration','FontSize',18);
    subplot(1,2,2);
    plot3(movingTest(:,1),movingTest(:,2),movingTest(:,3),'.r','MarkerSize',4); hold on;
    plot3(allCTControlPoints_FOV(:,1),allCTControlPoints_FOV(:,2),allCTControlPoints_FOV(:,3),'.b','MarkerSize',4); daspect([1 1 1]); grid on;
    title('After Fine Registration','FontSize',18);
    
%% save out the file
    transform_RAS = [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1]
    registrationMatrix_RAS = transform_RAS*registrationMatrix*transform_RAS;
    
    save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\combined_transformation_ini_LPS.txt', 'registrationMatrix', '-ascii');
    save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\combined_transformation_ini_RAS.txt', 'registrationMatrix_RAS', '-ascii');
    
    save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\vessel_transformation_ini_LPS.txt', 'vessel_transformationCPD_LPS', '-ascii');
    save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\PT-12\US\registration\centerline_transformationCPD_ini_LPS.txt', 'centerline_transformationCPD_LPS', '-ascii');
    

    
    %     save('E:\PROGRAM\Project_PhD\Registration\Data\MR_3DUS_healthy_study\CenterlineExtraction\LHV-05\MR_US_02\registration\transformation_ini_RAS.txt', 'registrationMatrix_RAS', '-ascii');





