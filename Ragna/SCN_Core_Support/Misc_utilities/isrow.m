% ISROW True if array is a ROW vector.
%	ISROW(V) returns logical true (1) if V is a 1 x n vector
%	where n >= 0, and logical false (0) otherwise.
%
%	See also iscol, isvector, isscalar, isnumeric
%		 islogical, ischar, isstrprop, isempty.

% created:
%	us	02-Feb-2006
% modified:
%	us	03-Feb-2006 23:46:41

function	tf=isrow(v)
		tf=isvector(v)&&size(v,1)==1;
