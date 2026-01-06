function y = vec(x)
%VEC Column-wise vectorization helper.
% Some toolboxes (and CVX internals) expect vec(x) to mean x(:).
% Keeping it here avoids path issues on systems where vec.m is missing.

y = x(:);
end
