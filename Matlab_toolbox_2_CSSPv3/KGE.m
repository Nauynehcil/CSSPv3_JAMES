function kge = KGE(sim, obs)
% KGE Calculates the Kling-Gupta Efficiency (KGE) metric
%   kge = KGE(sim, obs)
%
%   Input Requirements:
%       sim and obs must be vectors of the same length (row or column).
%
%   Output:
%       kge - Kling-Gupta Efficiency. A value closer to 1 indicates 
%             better model performance.
%
%   Note: This function depends on the custom 'CC.m' function.

% Check 1: Vector validation
if ~(isvector(sim) && isvector(obs))
    error('KGE:InputError', 'Inputs "sim" and "obs" must be vectors (n,1) or (1,n).');
end

% Check 2: Dimensional consistency
if ~isequal(size(sim), size(obs))
    if isequal(length(sim), length(obs))
        obs = reshape(obs, size(sim));  % Match obs shape to sim
    else
        error('KGE:SizeMismatch', 'The lengths of "sim" and "obs" do not match.');
    end
end

% Handle Missing Data (NaNs)
% Synchronize NaNs to ensure calculation consistency
sim(isnan(obs)) = nan;
obs(isnan(sim)) = nan;

% Calculate Means (ignoring NaNs)
ave_obs = mean(obs, 'omitnan');
ave_sim = mean(sim, 'omitnan');

% Calculate KGE Components:
% 1. Correlation (cc) - using custom CC function
cc = CC(sim, obs);

% 2. Bias Ratio (br) - ratio of the means
br = ave_sim / ave_obs;

% 3. Variability Ratio (rv) - ratio of the coefficients of variation
rv = (std(sim, 'omitnan') / ave_sim) / (std(obs, 'omitnan') / ave_obs);

% Calculate Final KGE Score
% Formula: KGE = 1 - sqrt( (cc-1)^2 + (br-1)^2 + (rv-1)^2 )
kge = 1 - sqrt((cc - 1)^2 + (br - 1)^2 + (rv - 1)^2);

end