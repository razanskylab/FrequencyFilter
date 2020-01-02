function sigMat = filtfiltj(b, a, sigMat)
  % same as filtfilt but accepts single precision inputs and only works on
  % 3D matrices!

  % Input validation
  % Only double precision is supported
  validateattributes(b, {'double'}, {'2d', 'finite', 'nonempty'}, 'filtfilt');
  validateattributes(a, {'double'}, {'2d', 'finite', 'vector', 'nonempty'}, 'filtfilt');

  sz = size(sigMat);

  % Parse SOS matrix or coefficients vectors and determine initial conditions
  [b, a, zi, nfact, L] = getCoeffsAndInitialConditions(b, a, sz(1));

  sigMat = reshape(sigMat, sz(1), []); %converting N-D matrix to 2-D.

  for ii = 1:L
    xt = bsxfun(@minus, 2 * sigMat(1, :), sigMat(nfact(1, 1) + 1:-1:2, :));
    [~, zo] = filter(b(:, ii), a(:, ii), xt, zi(:, ii) * xt(1, :)); % outer product
    [yc2, zo] = filter(b(:, ii), a(:, ii), sigMat, zo);
    xt = bsxfun(@minus, 2 * sigMat(end, :), sigMat(end - 1:-1:end - nfact(1, 1), :));
    yc3 = filter(b(:, ii), a(:, ii), xt, zo);
    [~, zo] = filter(b(:, ii), a(:, ii), yc3(end:-1:1, :), zi(:, ii) * yc3(end, :)); % outer product
    yc5 = filter(b(:, ii), a(:, ii), yc2(end:-1:1, :), zo);
    sigMat = yc5(end:-1:1, :);
  end

  sigMat = reshape(sigMat, sz); % To match the input size.

end

%--------------------------------------------------------------------------
function [b1, a1, zi, nfact, L] = getCoeffsAndInitialConditions(b, a, Npts)

  [L, ncols] = size(b);
  na = numel(a);

  % Rules for the first two inputs to represent an SOS filter:
  % b is an Lx6 matrix with L>1 or,
  % b is a 1x6 vector, its 4th element is equal to 1 and a has less than 2
  % elements.
  if ncols == 6 && L == 1 && na <= 2

    if b(4) == 1
      coder.internal.warning('signal:filtfilt:ParseSOS', 'SOS', 'G');
    else
      coder.internal.warning('signal:filtfilt:ParseB', 'a01', 'SOS');
    end

  end

  %----------------------------------------------------------------------
  % b and a are vectors that define the transfer function of the filter
  %----------------------------------------------------------------------

  coder.internal.errorIf((~isvector(a) ||~isvector(b)), 'signal:filtfilt:InputNotSupported');
  L = 1;
  % Check coefficients
  b1 = b(:);
  a1 = a(:);
  nb = numel(b);
  nfilt = max(nb, na);
  nfact = max(1, 3 * (nfilt - 1)); % length of edge transients

  % input data too short
  coder.internal.errorIf(Npts <= nfact(1, 1), 'signal:filtfilt:InvalidDimensionsDataShortForFiltOrder', nfact(1, 1));

  % Zero pad shorter coefficient vector as needed
  if nb < nfilt
    b1 = [b1; zeros(nfilt - nb, 1)];
  elseif na < nfilt
    a1 = [a1; zeros(nfilt - na, 1)];
  end

  % Compute initial conditions to remove DC offset at beginning and end of
  % filtered sequence.  Use sparse matrix to solve linear system for initial
  % conditions zi, which is the vector of states for the filter b(z)/a(z) in
  % the state-space formulation of the filter.
  if nfilt > 1
    rows = [1:nfilt - 1, 2:nfilt - 1, 1:nfilt - 2];
    cols = [ones(1, nfilt - 1), 2:nfilt - 1, 2:nfilt - 1];
    vals = [1 + a1(2, 1), a1(3:nfilt, 1).', ones(1, nfilt - 2), -ones(1, nfilt - 2)];
    rhs = b1(2:nfilt, 1) - b1(1, 1) * a1(2:nfilt, 1);
    zi = sparse(rows, cols, vals) \ rhs;
    % The non-sparse solution to zi may be computed using:
    %      zi = ( eye(nfilt-1) - [-a(2:nfilt), [eye(nfilt-2); ...
    %                                           zeros(1,nfilt-2)]] ) \ ...
    %          ( b(2:nfilt) - b(1)*a(2:nfilt) );
  else
    zi = zeros(0, 1);
  end

  zi = single(zi);

end
