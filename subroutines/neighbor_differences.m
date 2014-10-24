%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @name neighbor_differences.m
% @author Parker Mills, Ahrens Lab, Carnegie Mellon 2008
%
% @brief For a 3x3 Neighborhood, find difference between neighbors
% 
% ==================== INPUT PARAMETERS ======================
% @param    x           (2D Float) A 3x3 matrix
%
% ==================== RETURNED DATA =========================
% @return   difference  (Float)    Desired report value
%
% ==================== ASSUMPTIONS ======================
% @assume   None
%
% ==================== DISPLAYED PRODUCTS ====================
% @product  None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function difference = neighbor_differences(x)


%% Initialize
center = x(2,2);



%% Preferences
% Setting to be a 4-way maximum was found to be optimal for injections of SPIO in brain tissue
% Other settings may be better for different datasets
ways = 4; % Set to 4 or 8-pixel neighborhood
method = 'maximum'; %What to report: 'maximum'=Max of neighbor differences 'sum'=Sum of neighbor differences



%% Calculate 4-way neighbor differences
a = abs(x(1,2) - center);
b = abs(x(2,1) - center);
c = abs(x(2,3) - center);
d = abs(x(3,2) - center);
% Correct for overestimation of distance that comes from not accounting for fact that on unit circle, 0 == 2*Pi
if (a > pi) a = 2*pi-a; end
if (b > pi) b = 2*pi-b; end
if (c > pi) c = 2*pi-c; end
if (d > pi) d = 2*pi-d; end



%% Calculate 8-way neighbor differences, if requested in preferences
if(ways==8)
    e = abs(x(1,1) - center);
    f = abs(x(1,3) - center);
    g = abs(x(3,1) - center);
    h = abs(x(3,3) - center);
    % Correct for overestimation of distance that comes from not accounting for fact that on unit circle, 0 == 2*Pi
    if (e > pi) e = 2*pi-e; end
    if (f > pi) f = 2*pi-f; end
    if (g > pi) g = 2*pi-g; end
    if (h > pi) h = 2*pi-h; end
end



%% Report maximum or sum
if(strcmp(method,'maximum'))
    if(ways==4)
        temp=sort([a b c d]);
        difference=temp(4);
    end
    if(ways==8)
        temp=sort([a b c d e f g h]);
        difference=temp(8);
    end
end

if(strcmp(method,'sum'))
    if(ways==4) difference = a+b+c+d; end
    if(ways==8) difference = a+b+c+d+e+f+g+h; end
end



%%%%%%%%%%%%%
%   EOF
%%%%%%%%%%%%%