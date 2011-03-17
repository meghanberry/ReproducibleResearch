function print_matrix(x,varargin)
% print_matrix(x,col names cell array, row names cell)
%
% tor wager
% prints matrix values as tab delimited, 2 decimal places
%
% Example:
% t = [1 2; 3 4; 5 6];
% print_matrix(t,{'col1' 'col2'},{'row1' 'row2' 'row3');

% set up
s = size(x,2);
str = [repmat('%3.4f\t',1,s) '\n'];

if length(varargin) > 0
    colnames = varargin{1};
end

if length(varargin) > 1
    rownames = varargin{2};
    colnames = [{' '} colnames];
end 

% Print names
if length(varargin) > 0, 
    for i = 1:length(colnames)
        fprintf(1,'%s\t',colnames{i}), 
    end
    fprintf(1,'\n')
end

% print matrix
if length(varargin) > 1
    for i = 1:size(x,1)
        fprintf(1,'%s\t',rownames{i}),
        fprintf(1,str,x(i,:));
    end
else
    disp(sprintf(str,x'))
end




    


return
