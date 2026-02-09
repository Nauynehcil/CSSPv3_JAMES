function out_bias = BIAS(sim, obs)
% BIAS Calculates the Mean Bias
%   out_bias = BIAS(sim, obs)
%
%   Requirements: sim and obs must be vectors of the same dimension.
%   Output: Mean Bias.

% Check 1: Vector validation
if ~(isvector(sim) && isvector(obs))
    error('BIAS:InputError', 'sim and obs must be vectors of dimension (n,1) or (1,n)');
end

% Check 2: Dimension unification
if ~isequal(size(sim), size(obs))
    if isequal(length(sim), length(obs))
        obs = reshape(obs, size(sim));
    else
        error('BIAS:SizeMismatch', 'The lengths of sim and obs do not match; comparison is not possible');
    end
end

% Remove NaN values
obs(isnan(sim)) = nan;
sim(isnan(obs)) = nan;

% Bias calculation
out_bias = mean(sim - obs, 'omitnan');
end