%%
% load('E:\PROGRAM\Project_PhD\Registration\bcpd\bcpd-master\demo\bcpd-rigid\HV05normY.txt')
% load('E:\PROGRAM\Project_PhD\Registration\bcpd\bcpd-master\demo\bcpd-rigid\HV05y.txt')
% load('E:\PROGRAM\Project_PhD\Registration\bcpd\bcpd-master\demo\bcpd-rigid\HV05normX.txt')
% load('E:\PROGRAM\Project_PhD\Registration\bcpd\bcpd-master\demo\bcpd-rigid\HV05R.txt')
% load('E:\PROGRAM\Project_PhD\Registration\bcpd\bcpd-master\demo\bcpd-rigid\HV05t.txt')
% load('E:\PROGRAM\Project_PhD\Registration\bcpd\bcpd-master\demo\bcpd-rigid\HV05s.txt')




function [output, normlizationMatrix, mu, scale]= normalization(input, flagNormalizedbyOther, input2)
    %% normalize
    if flagNormalizedbyOther == 0
        
        % compute mean value: mu
        n=size(input,1);
        mu = sum(input)/(n);
        % compute scale: scale
        value = sum(sum((input  - mu).*(input - mu)));
        scale = sqrt(value/(n*3));
  
    else
         % compute mean value: mu
        n=size(input2,1);
        mu = sum(input2)/(n);
        % compute scale: scale
        value = sum(sum((input2  - mu).*(input2 - mu)));
        scale = sqrt(value/(n*3));

    end
         %normalization
        x_translate = input - mu;
        output = x_translate / scale;
        normlizationMatrix = [1/scale 0 0 -mu(1)/scale; 0 1/scale 0 -mu(2)/scale; 0 0 1/scale -mu(3)/scale; 0 0 0 1];
        
end




