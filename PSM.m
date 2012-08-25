%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file PSM.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon
% @brief Creates phase slope magnitude (PSM), angle (PSA), Angle Coherence (AC)
%        and Composite images based on provided phase image. 2nd and 3rd derivatives
%        were not found to be useful for analysis of phase changes in MR images, but
%        perhaps trying it with other image types or biological samples my prove of use.
%
% ==================== INPUT PARAMETERS ======================
% @param    image             (2D/3D Float)    Phase-unwrapped, optionally high-pass filtered image
% @param    dimensionality    (String)         '2d': Calculation performed slicewise in 2D 
%                                              '3d': Calculation performed in full 3D
%
% ==================== RETURNED DATA =========================
%
% ====Conventional PSM Products (Magnitude, x-,y-,z-components) ====
% @return  PSMrun.PSM.mag        (2D/3D Float)  Phase slope magnitude image   (*Main product* - Magnitude of phase change over space)
%          PSMrun.PSM.x          (2D/3D Float)  X-direction phase slope magnitude image   (X-component of phase change over space (1st derivative))
%          PSMrun.PSM.y          (2D/3D Float)  Y-direction phase slope magnitude image   (Y-component of phase change over space (1st derivative))
%          PSMrun.PSM.z          (2D/3D Float)  Z-direction phase slope magnitude image   (Z-component of phase change over space (1st derivative))
%
% ====Unpublished/Experimental PSM Products (Phase slope angle, angle coherence, and composite images)====
% @return  PSMrun.PSA        (2D/3D Float)  Phase slope angle image       (The angle (in radians) of the direction in which phase is changing over space)
% @return  PSMrun.AC         (2D/3D Float)  Phase angle coherence image   (Measure of how consistent the phase slope angle (PSA) is in a 3x3 pixel region. Calculated only for x-y plane, and not for z.)
% @return  PSMrun.Composite  (2D/3D Float)  Composite image               (Product of PSM.mag and AC images (PSM .* AC). In theory, highlights rapid spatial changes in phase that are all happening in the same direction!)
%
%
% ==================== ASSUMPTIONS ======================
% @assume Input data is phase data that has been unwrapped ([0, 2*pi] boundaries removed), and optionally high-pass filtered
%
% ==================== DEPENDENCIES ==================
% @depend myfun.m
% @depend
%
% ==================== EXAMPLE USAGE ========================
% PSMrun = PSM(image, dimensionality)
% PSMrun = PSM(unwrapped_phase, '2d')
% PSMrun = PSM(unwrapped_phase, '3d')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PSMrun = PSM(image, dimensionality)


%% Preferences
pref_diff_level = 1; %Which derivative to calculate (1 = first derivative, 2 = 2nd derivative, etc.)
pref_calculate_experimental = 0; % Calculating experimental images computationally intensive for large 3D volumes



%% Gather information and perform sanity checks

% Get dimensions
[x_dim y_dim z_dim] = size(image);

% If 2D image is given and user requested 3D PSM
if(ndims(image) == 2 && strcmp(dimensionality,'3d'))
    error('PSM: 3D computation requested for 2D image. Quitting.');
end

% Ensure dimensionality is either 2D or 3D
if(~strcmp(dimensionality,'2d') && ~strcmp(dimensionality,'3d'))
    error('PSM: dimensionality must be "2d" or "3d". Quitting.');
end



%% Compute difference (first derivative by default preferences) between neighboring voxels in 2D or 3D
% Also reduces floating point precision from double to single
PSMrun.PSM.x = single(diff(image, pref_diff_level, 1));
PSMrun.PSM.y = single(diff(image, pref_diff_level, 2));
if(strcmp(dimensionality,'3d'))
    PSMrun.PSM.z = single(diff(image, pref_diff_level, 3));
end



%% Preallocate results
if(strcmp(dimensionality,'2d'))
    PSMrun.PSM.mag   = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim);
    PSMrun.PSA       = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim);
    PSMrun.AC        = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim);
    PSMrun.Composite = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim);
else
    PSMrun.PSM.mag   = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim - pref_diff_level);
    PSMrun.PSA       = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim - pref_diff_level);
    PSMrun.AC        = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim - pref_diff_level);
    PSMrun.Composite = ones(x_dim-pref_diff_level, y_dim-pref_diff_level, z_dim - pref_diff_level);
end



%% Compute magnitude and angle for each row, filling [mag, angle]
for j1 = 1:x_dim - pref_diff_level
    for j2 = 1:y_dim - pref_diff_level
        
        if(strcmp(dimensionality,'3d'))
            for j3 = 1:z_dim - pref_diff_level
                PSMrun.PSM.mag(j1,j2,j3) = sqrt( PSMrun.PSM.x(j1,j2,j3).^2 + PSMrun.PSM.y(j1,j2,j3).^2 + PSMrun.PSM.z(j1,j2,j3).^2 );
            end
            PSMrun.PSA(j1,j2,:) = atan2(PSMrun.PSM.x(j1,j2,1:z_dim-pref_diff_level), PSMrun.PSM.y(j1,j2,1:z_dim-pref_diff_level)); % Angle is only done for X-Y plane and doesn't factor in Z
        else
            PSMrun.PSM.mag(j1,j2,:) = sqrt( PSMrun.PSM.x(j1,j2,1:z_dim).^2 + PSMrun.PSM.y(j1,j2,1:z_dim).^2);
            PSMrun.PSA(j1,j2,:) = atan2(PSMrun.PSM.x(j1,j2,1:z_dim), PSMrun.PSM.y(j1,j2,1:z_dim)); % Angle is only done for X-Y plane and doesn't factor in Z
        end
        
    end
end



%% Compute Angle Coherence (AC) and Composite images
if(pref_calculate_experimental)
    if(strcmp(dimensionality,'3d'))
        for j3 = 1:z_dim - pref_diff_level
            PSMrun.AC(:,:,j3) = -nlfilter(PSMrun.PSA(:,:,j3), [3 3], @neighbor_differences);
        end
    else
        for j3 = 1:z_dim
            PSMrun.AC(:,:,j3) = -nlfilter(PSMrun.PSA(:,:,j3), [3 3], @neighbor_differences);
        end
    end
    PSMrun.Composite = PSMrun.PSM.mag .* PSMrun.AC;
end



%%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%%%%%
