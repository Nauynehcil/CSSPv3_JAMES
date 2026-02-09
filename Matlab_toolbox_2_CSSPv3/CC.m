function R = CC(sim, obs)
% CC Calculates the Pearson Correlation Coefficient
%   R = CC(sim, obs)
%
%   Input Requirements:
%       sim and obs must be vectors of the same length.
%       Returns the correlation coefficient R.

% Check 1: Vector validation
if ~(isvector(sim) && isvector(obs))
    error('CC:InputError', 'Inputs "sim" and "obs" must be vectors (n,1) or (1,n).');
end

% Check 2: Dimensional consistency
if ~isequal(size(sim), size(obs))
    if isequal(length(sim), length(obs))
        obs = reshape(obs, size(sim));
    else
        error('CC:SizeMismatch', 'The lengths of "sim" and "obs" do not match.');
    end
end

% Handle Missing Data (NaNs)
% Ensure that a NaN in one vector results in a NaN in the same position of the other
sim(isnan(obs)) = nan;
obs(isnan(sim)) = nan;

% De-meaning (Center the data)
ccr1 = sim - mean(sim, 'omitnan');
ccr2 = obs - mean(obs, 'omitnan');

% Calculation using the centered data formula
% R = Σ( (sim-μ_s) * (obs-μ_o) ) / ( √(Σ(sim-μ_s)^2) * √(Σ(obs-μ_o)^2) )
R = sum(ccr1 .* ccr2, 'omitnan') / ...
    (sqrt(sum(ccr1.^2, 'omitnan')) * sqrt(sum(ccr2.^2, 'omitnan')));
end