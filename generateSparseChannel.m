
function [activeTaps, gains, h] = generateSparseChannel(M, K, minAbs, maxAbs)
%GENERATESPARSECHANNEL Generate K random nonzero channel coefficients.
%
% [activeTaps, gains, h] = generateSparseChannel(M, K, minAbs, maxAbs)
%
% M         : Total channel length
% K         : Number of active taps
% minAbs    : Minimum coefficient magnitude
% maxAbs    : Maximum coefficient magnitude
% activeTaps: Random active tap indices
% gains     : Random signed coefficients at those indices
% h         : Full M x 1 channel vector

    validateattributes(M, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(K, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative', '<=', M});

    if minAbs < 0 || maxAbs < minAbs
        error('Require 0 <= minAbs <= maxAbs.');
    end

    activeTaps = sort(randperm(M, K));

    magnitudes = minAbs + (maxAbs - minAbs) * rand(K, 1);
    signs = 2 * (rand(K, 1) >= 0.5) - 1;

    gains = magnitudes .* signs;

    h = zeros(M, 1);
    h(activeTaps) = gains;
end
    