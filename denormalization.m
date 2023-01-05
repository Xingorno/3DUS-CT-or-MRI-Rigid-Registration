function output = denormalization(input, flagNormalizedbyOther, input2)
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
        %% denormalization
        x_resize_inv = input*scale;
        output = x_resize_inv + mu;
        
end