% ISCOL True if array is a COLUMN vector.
%	ISCOL(V) returns logical true (1) if V is a n x 1 vector
%	where n >= 0, and logical false (0) otherwise.
%
%	See also isrow, isvector, isscalar, isnumeric
%		 islogical, ischar, isstrprop, isempty.

% created:
%	us	02-Feb-2006
% modified:
%	us	03-Feb-2006 23:46:44

function	tf=iscol(v)
		tf=isvector(v)&&size(v,2)==1;
