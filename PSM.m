%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file PSM.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon
%
% @brief Creates phase slope magnitude (PSM), phase slope angle (PSA), 
%        phase slope Angle Coherence (AC), and phase slope Composite images
%        based on provided phase image.
%
%        2nd and 3rd derivatives of phase slope were not found to be useful for analysis 
%        of phase changes in MR images, but perhaps trying it with other image types 
%        or biological samples may prove useful.
%
% ==================== INPUT PARAMETERS ======================
% @param    image             (2D/3D Float)    Phase image that has been phase-unwrapped and (optionally) high-pass filtered
%
% @param    dimensionality    (String)         '2d': Calculation performed slicewise in 2D 
%                                              '3d': Calculation performed in full 3D
%
% ==================== RETURNED DATA =========================
%
% ==== Published Products (phase slope magnitude and its x-, y-, and z-components) ====
% @return  PSMrun.PSM.mag   (2D/3D Float)  Phase slope magnitude image   (*Main product* - Magnitude of phase change over space)
%          PSMrun.PSM.x     (2D/3D Float)  X-direction phase slope magnitude image   (X-component of phase change over space (1st derivative))
%          PSMrun.PSM.y     (2D/3D Float)  Y-direction phase slope magnitude image   (Y-component of phase change over space (1st derivative))
%          PSMrun.PSM.z     (2D/3D Float)  Z-direction phase slope magnitude image   (Z-component of phase change over space (1st derivative))
%
% ==== Unpublished and Experimental Products (Phase slope angle, phase slope angle coherence, and phase slope composite images)====
% @return  PSMrun.PSA        (2D/3D Float)  Phase slope angle image       (The angle (in radians) of the direction in which phase is changing over space)
% @return  PSMrun.AC         (2D/3D Float)  Phase angle coherence image   (Measure of how consistent the phase slope angle (PSA) is in a 3x3 pixel region. Calculated only for x-y plane, and not for z.)
% @return  PSMrun.Composite  (2D/3D Float)  Composite image               (Product of PSM.mag and AC images (PSM .* AC). In theory, highlights rapid spatial changes in phase that are all happening in the same direction!)
%
%
% ==================== ASSUMPTIONS ======================
% @assume Input image data is a phase image that has undergone phase unwrapping (all [0, 2*pi] boundaries removed), and (optionally) high-pass filtered
%
% ==================== DEPENDENCIES ==================
% @depend None
%
% ==================== EXAMPLE USAGE ========================
% PSMrun = PSM(phantom_unwrapped_phase, '2d')
% PSMrun = PSM(brain_unwrapped_phase, '3d')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PSMrun = PSM(image, dimensionality)


%% Preferences %%%

% Which derivative to calculate (1 = first derivative, 2 = 2nd derivative, etc.)
pref_diff_level = 1; 

% Calculating the experimental and unpublished products described above is
% computationally intensive for large 3D volumes. Change this value to '0' 
% to skip generating these products.
pref_calculate_experimental = 1; 




%% Gather information and perform sanity checks %%%

% Get dimensions
[x_dim y_dim z_dim] = size(image);

% If 1D image is given
if(ndims(image) == 1)
    error('PSM: 1D images not supported. Quitting.');
end

% If 2D image is given and user requested 3D PSM
if(ndims(image) == 2 && strcmp(dimensionality, '3d'))
    error('PSM: 3D computation requested for 2D image. Quitting.');
end

% Ensure dimensionality is either 2D or 3D
if(~strcmp(dimensionality, '2d') && ~strcmp(dimensionality, '3d'))
    error('PSM: Dimensionality must be "2d" or "3d". Quitting.');
end




%% Compute difference (first derivative by default preferences) between neighboring voxels in 2D or 3D
% Also reduces floating point precision from double to single
PSMrun.PSM.x = single(diff(image, pref_diff_level, 1));
PSMrun.PSM.y = single(diff(image, pref_diff_level, 2));
if(strcmp(dimensionality, '3d'))
    PSMrun.PSM.z = single(diff(image, pref_diff_level, 3));
end



%% Preallocate result matrices
new_x_dim = x_dim - pref_diff_level;
new_y_dim = y_dim - pref_diff_level;
new_z_dim = z_dim - pref_diff_level;
if(strcmp(dimensionality,'2d'))
    PSMrun.PSM.mag   = ones(new_x_dim, new_y_dim, z_dim);
    PSMrun.PSA       = ones(new_x_dim, new_y_dim, z_dim);
    PSMrun.AC        = ones(new_x_dim, new_y_dim, z_dim);
    PSMrun.Composite = ones(new_x_dim, new_y_dim, z_dim);
else
    PSMrun.PSM.mag   = ones(new_x_dim, new_y_dim, new_z_dim);
    PSMrun.PSA       = ones(new_x_dim, new_y_dim, new_z_dim);
    PSMrun.AC        = ones(new_x_dim, new_y_dim, new_z_dim);
    PSMrun.Composite = ones(new_x_dim, new_y_dim, new_z_dim);
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



%%%%%%%%%%%%%%%%%
%    EOF
%%%%%%%%%%%%%%%%%
