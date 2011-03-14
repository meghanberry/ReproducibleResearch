% ISMATRIX: Returns 1 if the input matrix is 2+ dimensional, 0 if it is a scalar 
%           or vector.
%
%     Usage ismat = ismatrix(X)
%

% RE Strauss, 5/19/00

function ismat = ismatrix(X)
  [r,c] = size(X);
  if (r>1 & c>1)
    ismat = 1;
  else
    ismat = 0;
  end;

  return;
