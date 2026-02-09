function mean_data = weightmean(data, area)
% WEIGHTMEAN Performs area-weighted averaging for 2D or 3D datasets
%   mean_data = weightmean(data, area)
%
%   Input:
%       data - Numeric array of 2D (lat x lon / lon x lat) 
%              or 3D (lat x lon x time / lon x lat x time).
%              Note: Pre-masking is required. Set non-weighted regions to NaN
%       area - Grid area matrix. Must be a 2D matrix matching the 
%              first two dimensions of "data".
%
%   Output:
%       mean_data - Scalar (if data is 2D) or a time-series vector 
%                   (if data is 3D).

% Check input dimensions
% Ensures "area" is 2D and matches the spatial extent of "data"
if ~ismatrix(area) || ~isequal(size(data, 1:2), size(area))
    error('Input "area" must be a 2D matrix matching the first two dimensions of "data".');
end

% Calculate total effective area for normalization
% Identifies valid (non-NaN) cells in the dataset
valid = ~isnan(data); 
total_area = sum(area .* valid, [1,2], 'omitnan');

% Prevent division by zero if no valid area is found
total_area(total_area == 0) = NaN;

% Weighted averaging calculation
if ndims(data) == 2
    % For 2D spatial data
    weight = area ./ total_area;
    mean_data = sum(data .* weight, 'all', 'omitnan');
elseif ndims(data) == 3
    % For 3D spatio-temporal data
    % Weight is normalized by total valid area per time step
    weight = area ./ total_area;
    mean_data = sum(data .* weight, [1,2], 'omitnan');
    mean_data = squeeze(mean_data); % Remove singleton dimensions
else
    error('Input "data" must be a 2D or 3D array.');
end

end