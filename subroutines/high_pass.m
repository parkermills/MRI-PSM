%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name high_pass.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon
% @brief High-pass filters an image - more frequently used than high_pass2
% 
% ==================== INPUT PARAMETERS ======================
% @param    image    (2D/3D Float)    Image to be high-pass filtered
% @param    size     (Float)          Filter size, in pixels, for frequency sampling
% @param    min      (Float)          Percentage of frequencies to be filtered out - range: [0.0, 1.0]
%
% ==================== RETURNED DATA =========================
% @return   result   (2D/3D Float)    High-passed result image
%
% ==================== ASSUMPTIONS ======================
% @assume   
%
% ==================== EXAMPLE USAGE =========================
% differences = calc_diff(im1,im2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result=high_pass(image,filter_size,min)


%% Find out if image is 2D or 3D
x = ndims(image);
if(x == 3)
    [x_dim y_dim z_dim] = size(image);
end


%% Create 2D filter
[f1,f2] = freqspace(filter_size,'meshgrid');
r = sqrt(f1.^2 + f2.^2);
Hd = ones(filter_size);
Hd(r < min) = 0;
h = fsamp2(Hd);


%% Apply high-pass filter to volume
% If image is 2D
if(x == 2)
    result = filter2(h,image,'same');
else
    % If image is 3D, high_pass each slice in same fashion
    if(x == 3)
        for j_1 = 1:z_dim
            result(:,:,j_1)=filter2(h,image(:,:,j_1),'same');
        end
    else
        disp('high_pass: Image must be 2D or 3D. Quitting.');
        return
    end
end

%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%