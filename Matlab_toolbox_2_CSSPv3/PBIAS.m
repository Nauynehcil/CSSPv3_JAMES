function out_pbias = PBIAS(sim, obs)
% PBIAS calculates the Percent Bias
%   out_pbias = PBIAS(sim, obs)
%
%   Requirements:
%       sim and obs must be vectors of the same dimension (n,1) or (1,n)
%
%   Output:
%       out_pbias - Percent Bias (unit: %)
%                   Positive values indicate overestimation, negative values indicate underestimation

% Check 1: Verify inputs are vectors
if ~(isvector(sim) && isvector(obs))
    error('PBIAS:InputError', 'sim and obs must be vectors of dimension (n,1) or (1,n)');
end

% Check 2: Unify dimensions
if ~isequal(size(sim), size(obs))
    if isequal(length(sim), length(obs))
        obs = reshape(obs, size(sim));
    else
        error('PBIAS:SizeMismatch', 'The lengths of sim and obs do not match');
    end
end

% Remove NaNs
sim(isnan(obs)) = nan;
obs(isnan(sim)) = nan;

% PBIAS Calculation
out_pbias = 100 * sum(sim - obs, 'omitnan') / sum(obs, 'omitnan');
end