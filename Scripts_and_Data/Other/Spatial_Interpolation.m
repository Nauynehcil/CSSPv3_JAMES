%% MATLAB Example: Spatial Interpolation and NetCDF Export
clc; clear;

% Define study period
st_year = 1991;
en_year = 2020;
tl = (en_year - st_year + 1) * 12;

% Load source and target grid coordinates
% Read source Lon/Lat (Also as example of 1D coordinate vectors)
lon_src = ncread('source_data.nc', 'lon');
lat_src = ncread('source_data.nc', 'lat');
[xlat_src, xlon_src] = meshgrid(lat_src, lon_src);

% Read target Lon/Lat (Also as example of 2D coordinate matrices)
lon_target = ncread('target_grid.nc', 'lon');
lat_target = ncread('target_grid.nc', 'lat');

% Pre-allocate 3D memory buffers for processed data (Lon x Lat x Time)
% Optimized for memory efficiency during long-term time series processing
processed_ET  = zeros(length(lon_target), length(lat_target), tl);
processed_SMs = zeros(length(lon_target), length(lat_target), tl);

n = 1;
for year = st_year:en_year 
    % Load raw monthly data from annual NetCDF files
    et_raw  = ncread(['Origin_Data_', num2str(year), '.nc'], 'ET'); 
    sms_raw = ncread(['Origin_Data_', num2str(year), '.nc'], 'SMs');
    
    for mon = 1:12
        % Core Logic: Bilinear interpolation to unify spatial resolution
        % Maps source grid values to the target grid resolution using 'linear' method
        processed_ET(:,:,n)  = interp2(xlat_src, xlon_src, et_raw(:,:,mon), lat_target, lon_target, 'linear');
        processed_SMs(:,:,n) = interp2(xlat_src, xlon_src, sms_raw(:,:,mon), lat_target, lon_target, 'linear');
        n = n + 1;
    end
    fprintf('Year %d processed.\n', year);
end

% 3. Export consolidated dataset to NetCDF4 format
output_file = 'Consolidated_Dataset.nc';
ncid = netcdf.create(output_file, 'NETCDF4');

% Define NetCDF dimensions: Longitude, Latitude, and Time
dimX = netcdf.defDim(ncid, 'lon', length(lon_target));
dimY = netcdf.defDim(ncid, 'lat', length(lat_target));
dimT = netcdf.defDim(ncid, 'time', n - 1);

% Define variables and configure global FillValue (NaN)
var_names = {'ET', 'SMs'};
var_data  = {processed_ET, processed_SMs};

for i = 1:length(var_names)
    % Define variable metadata
    varid = netcdf.defVar(ncid, var_names{i}, 'double', [dimX dimY dimT]);
    
    % Set FillValue to ensure NaN handling consistency
    netcdf.defVarFill(ncid, varid, false, nan);
    
    % Switch to data mode for writing
    netcdf.endDef(ncid); 
    netcdf.putVar(ncid, varid, var_data{i});
    
    % Switch back to define mode for the next variable iteration
    if i < length(var_names)
        netcdf.redef(ncid);
    end
end

% Close file and release system handle
netcdf.close(ncid);
disp('NetCDF export completed successfully.');