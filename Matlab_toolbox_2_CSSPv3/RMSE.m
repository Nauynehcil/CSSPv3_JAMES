function out_rmse = RMSE(sim, obs)
% RMSE Calculates the Root Mean Square Error
%   out_rmse = RMSE(sim, obs)
%
%   Input Requirements:
%       sim and obs must be vectors of the same dimension (n,1) or (1,n).
%
%   Output:
%       out_rmse - The Root Mean Square Error between simulated and 
%                  observed values. Lower values indicate better performance.

% Check 1: Vector validation
if ~(isvector(sim) && isvector(obs))
    error('RMSE:InputError', 'Inputs "sim" and "obs" must be vectors (n,1) or (1,n).');
end

% Check 2: Dimensional consistency
if ~isequal(size(sim), size(obs))
    if isequal(length(sim), length(obs))
        obs = reshape(obs, size(sim));
    else
        error('RMSE:SizeMismatch', 'The lengths of "sim" and "obs" do not match.');
    end
end

% Handle Missing Data (NaNs)
% Synchronize NaNs between simulated and observed data
sim(isnan(obs)) = nan;
obs(isnan(sim)) = nan;

% Calculation of RMSE
% Formula: RMSE = sqrt( mean( (sim - obs)^2 ) )
out_rmse = sqrt(mean((sim - obs).^2, 'omitnan'));

end