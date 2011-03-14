function img = scn_map_image(loadImg, sampleTo, varargin)
    % [imgData] = scn_map_image(loadImg, sampleTo, varargin)
    %
    % Tor Wager, July 2007
    %
    % This function takes an image name in loadImg
    % and loads the data, resampling to the space defined
    % in the image sampleTo.
    %
    % Optional:
    % 'write' followed by name of resampled image to write
    %
    % Compatible with SPM5, don't know yet about other versions.
    % One image file in each.  Don't know how it'll handle filenames with
    % comma, volume number at end.
    %
    % Examples:
    % img = scn_map_image(EXPT.mask, EXPT.SNPM.P{1}(1,:), 'write', 'resliced_mask.img');
    %
    %

    % optional input arguments
    % -------------------------------------------------------
    if ~isempty(varargin)
        for i = 1:length(varargin)
            if ischar(varargin{i})
                switch varargin{i}

                    % functional commands
                    case {'write', 'name', 'outname'}
                        outname = varargin{i+1};
                        varargin{i+1} = [];

                    otherwise, warning(['Unknown input string option:' varargin{i}]);
                end
            end
        end
    end


    % image whose space to sample into
    % -------------------------------------------------------
    sampleTo = deblank(sampleTo);
    Vto = spm_vol(sampleTo);               % volume to sample TO

    Vto = Vto(1); % for SPM5 compat
    
    if isempty(Vto)
        fprintf('%s is empty or missing.\n', sampleTo)
        return
    end

    Mto = Vto(1).mat;                        % mat to sample TO


    % image to load and sample to new space
    % -------------------------------------------------------
    loadImg = deblank(loadImg);
    Vmap = spm_vol(loadImg);                % mask to resample
    if isempty(Vmap)
        fprintf('%s is empty or missing.\n', loadImg)
        return
    end

    Mmap = Vmap(1).mat;



    % image dimensions
    % -------------------------------------------------------
    dim = Vto(1).dim(1:3);


    % output image data
    % -------------------------------------------------------
    img = zeros(dim);


    % map slice-by-slice
    % -------------------------------------------------------
    for j = 1:dim(3)

        Mslice  = spm_matrix([0 0 j]);      % Matrix specifying this slice

        Mtrans  = Mto \ Mmap \ Mslice;          % Affine mappping mtx: Mask -> TOvol

        img(:, :, j) = spm_slice_vol(Vmap(1), Mtrans, dim(1:2), [0 NaN]);


    end

    % create image, if asked for
    % -------------------------------------------------------
    if exist('outname', 'var')
        Vout = Vto;
        Vout.fname = outname;
        Vout  = spm_create_vol(Vout);

        spm_write_vol(Vout, img);
    end



end
