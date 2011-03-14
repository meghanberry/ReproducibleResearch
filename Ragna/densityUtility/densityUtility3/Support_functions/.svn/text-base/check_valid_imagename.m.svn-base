function P = check_valid_imagename(P, varargin)
% P = check_valid_imagename(P, [return error])
%
% optional input: 1 to break with an error message, or nothing to select
% missing filenames graphically

for i = 1:size(P,1)

    img = deblank(P(i, :));

    % take off trailing commas (SPM notation for multiple vols)
    wh = find(img == ',');
    if ~isempty(wh)
        wh = wh(end);
        img(wh:end) = [];
    end

    if ~(exist(img, 'file'))

        if ~isempty(varargin) && varargin{1}
            error(sprintf('Cannot find image:\n%s\n', P(i, :)));


        else
            fprintf(1,'Cannot find image:\n%s\nPlease select correct directory.\n',P(i,:));

            dd = spm_select(1,'dir','Select directory.');
            dd = [dd filesep];

            [myd] = fileparts(P(i,:));
            len = length(myd);

            P = [repmat(dd,size(P,1),1) P(:,len+1:end)];

            disp(P(i,:));

        end

    end
end


return

