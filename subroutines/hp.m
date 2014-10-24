%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @file hp.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon
% @brief Quick wrapper for high-pass filtering with good parameters
%
% ==================== INPUT PARAMETERS ======================
% @param image         (2D/3D float)      Image
%
% ==================== RETURNED DATA =========================
% @return hped         (2D/3D float)      High-passed image
%
% ==================== ASSUMPTIONS ===========================
% @assume  
%
% ==================== EXAMPLE USAGE =========================
% hped = hp(im1,im2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hped = hp(image)

% Paramters found to be good
%(d = 5, frequency = 0.1)
% High-pass freq_cutoff, when set to 0.1, is designed to remove lowest 10% of frequencies.
% Setting freq_cutoff > 0.1 provides no extra benefit. Anything smaller has no effect.
% Filter sizes > 5 were used to no extra benefit. Filter sizes <=3 are too small, and produce artifacts.

% Preferences
pref_high_pass_filter_size = 5;   % High-pass filter size in pixels %used to be 5
pref_high_pass_freq_cutoff = 0.1; %0.06 % Value in range [0.0, 1.0] that determines % of frequencies to remove

% High pass operation
hped = high_pass(image, pref_high_pass_filter_size, pref_high_pass_freq_cutoff);

%%%%%%%%%%%%%%%
%     EOF
%%%%%%%%%%%%%%%