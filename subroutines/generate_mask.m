%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name generate_mask.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2009
% @brief Creates mask from MRI magnitude image by fitting intensity values
%        to either a Rayleigh-Gaussian distribution or a
%        Rayleigh-Gaussian-Gaussian (dual Gaussian) distribution.
%
% ==================== INPUT PARAMETERS ======================
% @param MRIdata         (3D Float)     MRIdata object
% @param mag_threshold   (Float)     Multiple of standard deviation. Examples for Rayleigh distribution include:
%                                       1.0 = 39% of noise values
%                                       2.0 = 86% of noise values
%                                       2.5 = 95.6% of noise values
%                                       3.0 = 98.8% of noise values
%
% @param phase_threshold (Float)     Multiple of standard deviation.
% @param dual_gaussian   (Float)        1.0 = yes, dual gaussian; 0.0 = no, single gaussian (default)
%
%
%==================== RETURNED DATA =========================
% @return mask          (3D Logical)    Conventional Image mask (Binary logical)
%
% ==================== ASSUMPTIONS ===============================
% @assume MRI magnitude noise has a Rician noise distribution
% @assume MRI magnitude signals have Gaussian signal distributions
% @assume Signal-to-Noise Raio (SNR) of dataset is greater than 2.0 for proper functioning
%
% ==================== SAMPLE USAGE ==========================
% mask = generate_mask(mag, phase_threshold, mag_threshold, dual_gaussian);
% mask = generate_mask(mag, phase_threshold, mag_threshold, dual_gaussian);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Rayleigh Fundamentals
% PDF = (x/sigma^2)   exp( -x^2 / (2 sigma^2) )
% mean = sigma * sqrt(pi/2)
% mode = sigma
% var = (4-pi) * sigma / 2
% [0      sigma] = 39% of values
% [0  2.0*sigma] = 86% of values
% [0  2.5*sigma] = 95.6% of values
% [0  3.0*sigma] = 98.8% of values

%% Gaussian Fundamentals
% PDF = 1/(sigma sqrt(2pi)) exp (- (x-u)^2 / (2 sigma^2) )
% mean = u
% mode = u
% var = sigma^2
% [-inf      sigma] = 84.2% of values
% [-inf  2.0*sigma] = 97.8% of values
% [-inf  3.0*sigma] = 99.9% of values

function mask = generate_mask(MRIdata, mag_threshold, phase_threshold, dual_gaussian, showfigures)


%% Preferences
fillholes = 1;
removeislands = 1;
pref_pixels_shrink_mask = 1.0; % Number of pixels to shrink mask when done filling holes
pref_bin_fraction = 0.01; % bins per pixel (e.g., 0.02 sets num_bins to be 2% of the number of pixels)
% For large datasets, set to 0.00005. Small datasets, 0.01
pref_detailed_fit_summary = 1; % Display a detailed fit summary?
pref_phase_masking = 0; % Should we also perform phase masking?


%% Initializations

% If this isn't an MRI data file
if(~isfield(MRIdata,'mag'))
    MRIdata.mag = MRIdata;
end

% Make complex data into magnitude data
if(~isreal(MRIdata.mag))
    MRIdata.mag = abs(MRIdata.mag);
end


% Choose between dual- or single Gaussian
if(~exist('dual_gaussian','var'))
    dual_gaussian = 0;
end

% Show user what they asked for
disp(['generate_mask: Chosen magnitude threshold (factor of sigma) is ',num2str(mag_threshold)]);
disp(['generate_mask: Chosen phase threshold (factor of sigma) is ',num2str(phase_threshold)]);

% Store dimensions and maximum image values
if(ndims(MRIdata.mag) == 3)
    [x_dim y_dim z_dim] = size(MRIdata.mag);
    max_mag = max(max(max(MRIdata.mag)));
end
if(ndims(MRIdata.mag) == 2)
    [x_dim y_dim] = size(MRIdata.mag);
    z_dim = 1;
    max_mag = max(max(MRIdata.mag));
end

% Calculate the number of pixels and histogram bins
pixels = x_dim * y_dim * z_dim;
num_bins = pixels * pref_bin_fraction;

% Create histogram to summarize intensity distribution
mag_hist = hist(abs(reshape(MRIdata.mag, 1 ,x_dim*y_dim*z_dim)),num_bins).';


%% Fit Initializations
% Order of parameters in vectors is: [ a1 a2 a3      b2 b3    c1 c2 c3   d ]
%                                    [ AMPLITUDEs   MEANs    STDEVs    OFFSET ]

% Set lower bounds
% Rationale: For magnitude images, all lower bounds are 0.0
lower_bound      = [0.0, 0.0,         0.0,         0.0, 0.0,         0.0];
lower_bound_dual = [0.0, 0.0, 0.0,    0.0, 0.0,    0.0, 0.0, 0.0,    0.0];

% Set upper bounds
% Rationale: noise_mean_amplitude & signal_mean_amplitude must be < pixels
%            signal_mean must be < num_bins
%            noise_stdev & signal_stdev must be < num_bins/2.0
%            const_offset must be < pixels/10.0
upper_bound      = [pixels, pixels,                 num_bins,               num_bins/2.0, num_bins/2.0,                   pixels/10.0];
upper_bound_dual = [pixels, pixels, pixels,         num_bins, num_bins,     num_bins/2.0, num_bins/2.0, num_bins/2.0,     pixels/10.0];

% Set initial values
% Rationale: noise_mean_amplitude & signal_mean_amplitude will be on the order of pixels
%            signal_mean is likely to be near the center of the distribution
%            noise_stdev & signal_stdev are likely to be around 1/30 the number of bins
%            constant_offset is likely very small, so I made it 10
initial_values      = [pixels  pixels/16.0                        (num_bins/4.0)                          (num_bins/30.0)  (num_bins/3.0)                      10];
initial_values_dual = [pixels  pixels/16.0  pixels/16.0           (num_bins/4.0)  (num_bins/2.0)          (num_bins/30.0)  (num_bins/3.0)   (num_bins/3.0)     10];

% Set fit options - these fits assume a Rayleigh noise distribution (mean = 0.0) and either one, or two, Gaussian signal distributions
fit_options      = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', lower_bound,      'Upper', upper_bound,      'Startpoint', initial_values     );
fit_options_dual = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', lower_bound_dual, 'Upper', upper_bound_dual, 'Startpoint', initial_values_dual);
fit_type      = fittype('a1 * x * exp(-(x^2)/(2*(c1^2))) / (c1^2)       +      a2 * (1 / (c2 * sqrt(2 * pi))) * exp(-((x-b2)^2)/(2*(c2^2)))                                                                          + d','options',fit_options);
fit_type_dual = fittype('a1 * x * exp(-(x^2)/(2*(c1^2))) / (c1^2)       +      a2 * (1 / (c2 * sqrt(2 * pi))) * exp(-((x-b2)^2)/(2*(c2^2)))      +      a3 * (1 / (c3 * sqrt(2 * pi))) * exp(-((x-b3)^2)/(2*(c3^2))) + d','options',fit_options_dual);
x_mag = [1:num_bins].';



%% Run fit and extract values
% Run fit
[fit_summary,    fit_quality,    fit_results  ] = fit(x_mag, mag_hist, fit_type  );
[fit_summary_dual,  fit_quality_dual,  fit_results_dual] = fit(x_mag, mag_hist, fit_type_dual);

% If preferred, display fit information
if(pref_detailed_fit_summary)
    disp(fit_summary); % Displays all detailed fit information
    disp(fit_summary_dual);
end

% Extract fit values
fit_summary_values = coeffvalues(fit_summary);
fit_summary_values_dual = coeffvalues(fit_summary_dual);

a1 = fit_summary_values(1);
a2 = fit_summary_values(2);
b2 = fit_summary_values(3);
c1 = fit_summary_values(4);
c2 = fit_summary_values(5);
d  = fit_summary_values(6);

a1_dual = fit_summary_values_dual(1);
a2_dual = fit_summary_values_dual(2);
a3_dual = fit_summary_values_dual(3);
b2_dual = fit_summary_values_dual(4);
b3_dual = fit_summary_values_dual(5);
c1_dual = fit_summary_values_dual(6);
c2_dual = fit_summary_values_dual(7);
c3_dual = fit_summary_values_dual(8);
d_dual  = fit_summary_values_dual(9);


%% Perform SNR calculations and display info

% Determine ordering of peak mean values for single-Gaussian fit
mean_rayleigh = c1 * sqrt(pi/2);
mean_gaussian = b2;
if(mean_rayleigh < mean_gaussian)
    noise_mean = mean_rayleigh;
    noise_sigma = c1;
    SNR = mean_gaussian / mean_rayleigh
else
    noise_mean = mean_gaussian;
    noise_sigma = c2;
    SNR = mean_rayleigh / mean_gaussian
end
hist_threshold = noise_mean + mag_threshold * noise_sigma; % Note: sigma means different things for these 2 different distributions


% Determine ordering of peak mean values for dual-Gaussian fit
% Ensure that central peak, and lowest peak, are the ones that are thresholded out!
% Not the main signal, and not just the lowest peak.
mean_rayleigh_dual = c1_dual * sqrt(pi/2);
mean_gaussian_1_dual = b2_dual;
mean_gaussian_2_dual = b3_dual;
sigmas = [c1_dual c2_dual c3_dual];
offsets = [0 mean_gaussian_1_dual mean_gaussian_2_dual];
[sorted sorted_ix] = sort([mean_rayleigh_dual mean_gaussian_1_dual mean_gaussian_2_dual]);
hist_threshold_dual = offsets(sorted_ix(2)) + mag_threshold * sigmas(sorted_ix(2));

% % Display SNR information
% disp(['generate_mask: SNR (signal_mean / noise_var) = ',num2str(SNR_var)]);
% disp(['generate_mask: CNR ((signal_mean-noise_mean) /noise_var) = ',num2str(CNR_var)]);
% disp(['generate_mask: SNR (signal_mean / noise_mean) = ',num2str(SNR_mean)]);
% disp(['generate_mask: CNR ((signal_mean-noise_mean)/noise_mean) = ',num2str(CNR_mean)]);
% disp(['generate_mask: "Quality" (signal_mean-noise_mean)/sqrt(signal_var + noise_var) = ',num2str((b2-noise_mean)/sqrt(noise_var + c2^2))]);



%% Calculate true threshold
if(dual_gaussian)
    true_threshold = (hist_threshold_dual / num_bins) * max_mag;
else
    true_threshold = (hist_threshold / num_bins) * max_mag;
end

%true_phase_threshold = phase_threshold/SNR_var;
disp(['generate_mask: True magnitude threshold is ', num2str(true_threshold)]);
%disp(['generate_mask: True phase threshold is ',     num2str(true_phase_threshold)]);



%% Plot magnitude image histogram, noise-signal fit, and masking threshold to allow user sanity-check

if(showfigures)
    figure;
    title('All fits. Blue = Data, Red = Fit, Green = Threshold Cutoff for Mask Generation');
    xlabel('Pixel Intensity');
    ylabel('Amplitude (relative pixel number)');
    hold on;
    
    % Plot magnitude in blue
    plot(mag_hist,'b','LineWidth',2);
    
    % Plot fit in red
    if(~dual_gaussian)
        plot(a1   .* x_mag .* exp(-(x_mag .^ 2)./(2.0 .* (c1   .^ 2))) / (c1  .^2) + a2   .* (1.0 / (c2   .* sqrt(2.0 * pi))) .* exp(-((x_mag-b2  ).^2)/(2.0 .* (c2  .^2)))                                                                                  + d,'r','LineWidth',2);
    else
        plot(a1_dual .* x_mag .* exp(-(x_mag .^ 2)./(2.0 .* (c1_dual .^ 2))) / (c1_dual.^2) + a2_dual .* (1.0 / (c2_dual .* sqrt(2.0 * pi))) .* exp(-((x_mag-b2_dual).^2)/(2.0 .* (c2_dual.^2))) + a3_dual .* (1.0 / (c3_dual .* sqrt(2.0 * pi))) .* exp(-((x_mag-b3_dual).^2)/(2.0 .* (c3_dual.^2))) + d_dual,'r','LineWidth',2);
    end
    
    % Plot mask threshold in green
    if(~dual_gaussian)
        plot(a1 ./ c1 ./ 2.0 .*  (x_mag < hist_threshold),'g','LineWidth',3); % A bright green line to show the chosen cutoff threshold
    else
        plot(a1 ./ c1 ./ 2.0 .*  (x_mag < hist_threshold_dual),'g','LineWidth',3); % A bright green line to show the chosen cutoff threshold
    end
    hold off;
    
end

%% Generate Mask using calculated threshold
mask = logical(MRIdata.mag > true_threshold);
images(mask);


%% Perform phase masking, if preferred (DISABLED IN PREFERENCES UNTIL LATER NOTICE)
if(pref_phase_masking)
    phase_mask = logical((MRIdata.phase > true_phase_threshold ) & (MRIdata.phase < -true_phase_threshold));
    combined_mask = logical(mask & phase_mask);
    
    % Go through combined mask, checking for connectivity
    Tm = 4;
    Tp = 4;
    for j_1=2:x_dim-1
        for j_2=2:y_dim-1
            
            if (z_dim > 2)
                for j_3=2:z_dim-1
                    if(~combined_mask(j_1,j_2,j_3))
                        
                        % if >Tm voxels border in mag mask, then restore pixel
                        if(sum(sum(mask(j_1-1:j_1+1,j_2-1:j_2+1,j_3-1:j_3+1))) > Tm)
                            combined_mask(j_1,j_2,j_3) = 1.0;
                        end
                        % if >Tp voxels border in phase mask, then restore pixel
                        if(sum(sum(phase_mask(j_1-1:j_1+1,j_2-1:j_2+1,j_3-1:j_3+1))) > Tp)
                            combined_mask(j_1,j_2,j_3) = 1.0;
                        end
                    end
                end
            end
            
            if (z_dim == 1)
                if(~combined_mask(j_1,j_2,1))
                    
                    % if >Tm voxels border in mag mask, then restore pixel
                    if(sum(sum(mask(j_1-1:j_1+1,j_2-1:j_2+1))) > Tm)
                        combined_mask(j_1,j_2,1) = 1.0;
                    end
                    % if >Tp voxels border in phase mask, then restore pixel
                    if(sum(sum(phase_mask(j_1-1:j_1+1,j_2-1:j_2+1))) > Tp)
                        combined_mask(j_1,j_2,1) = 1.0;
                    end
                end
            end
            
        end
    end
end


%% Zero out the border
mask(1,:) = 0;
mask(z_dim,:) = 0;
mask(:,1) = 0;
mask(:,z_dim) = 0;


%%
if(removeislands)
    maskshrunk = mask;
    for j_1=2:x_dim-1
        for j_2=2:y_dim-1
            for j_3=2:z_dim-1
                if(mask(j_1,j_2,j_3))
                    if(~(   mask(j_1+1,j_2,j_3) || mask(j_1-1,j_2,j_3) || mask(j_1,j_2+1,j_3) || mask(j_1,j_2-1,j_3) || mask(j_1,j_2,j_3+1) || mask(j_1,j_2,j_3-1)    ))
                        maskshrunk(j_1,j_2,j_3) = 0.0;
                    end
                end
            end
        end
    end
    mask = maskshrunk;
end



%% Fill holes present in mask, then region shrink by preferred # of pixels
if(fillholes)
    mask = imfill(mask,'holes'); % Fill holes
    maskshrunk = mask;
    for k_1=1:pref_pixels_shrink_mask % How many pixels you would like to shrink mask by
        for j_1=2:x_dim-1
            for j_2=2:y_dim-1
                for j_3=2:z_dim-1
                    if(~mask(j_1,j_2,j_3))
                        maskshrunk(j_1-1:j_1+1,j_2-1:j_2+1,j_3-1:j_3+1) = 0.0;
                    end
                end
            end
        end
        mask = maskshrunk;
    end
end



%% Show various images for sanity check
if(showfigures)
    images(MRIdata.mag); title(['generate_mask: Provided image']);
    images(mask); title(['generate_mask: Magnitude mask generated from provided image']);
end


% If we're performing phase masking as well
if(pref_phase_masking)
    images(phase_mask); title(['generate_mask: Phase mask generated from provided image']);
    images(combined_mask); title(['generate_mask: Combined mask, without connectivity']);
    images((combined_mask) .* MRIdata.mag); title(['generate_mask: Combined mask on mag image']);
end



%%%%%%%%%%%%%%%%%%%%%
% EOF
%%%%%%%%%%%%%%%%%%%%%