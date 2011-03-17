function [model,delta] = getPredictors(stimList, HRF,varargin)
% function [model,delta] = getPredictors(stimList, HRF, [downsample factor])
%
% Build predictors and delta functions, given a condition function or delta
% function and either a convolution matrix or vector.
%
% IMPORTANT: YOU MUST ADD THE INTERCEPT YOURSELF!
%
% stimList: condition function or delta function
% HRF:      1) hemodynamic response function
%           2) Basis set (columns)
%           3) or convolution matrix (columns
%               are HRF), defined as:
%               HRF = tril(toeplitz(hrf));
%
%           multiple column vectors for HRF are treated as basis functions!
% [downsamp]    takes every nth element of the design matrix
%
%inputs:
% 1     a col. vector of stimulus conditions OR a delta function matrix
% 2	    an HRF vector sampled at the frequency of the stimulus vector, OR
%	    a convolution matrix H
%outputs:
% 1     a n x 2 matrix of regressors (cols) for each condition
% 2     a n x k delta matrix with onsets
%
% Tor Wager, last modified 2/22/04 to center predictors
%            modified to optionally take H convolution matrix as input
%            in place of HRF vector.  See convmtx.m (Buracas,mseq toolbox)
%
% stimList can be condition function e.g., [1 3 2 4 3 2 1]' or
% delta matrix (n x k), n samples and k conditions, e.g., [1 0 0 0 1 0 1]'
%
% Resampling: the default N in matlab resample has built-in antialiasing,
% but may not be good for fmri designs!  The appropriate downsampling
% is expected to be res*TR (res is units of samples/s), but we use 0
% because the model will depend on the analysis method used, and this is
% the most veridical approach.  With N = 0, every ith sample is used, where
% i is the downsampling factor you input.  Popular choices are 16*TR (for
% onsets2delta.m), using the SPM default res of 16.
% Delta is NOT resampled.
%
% example: TR = 2, 16 samples per second in hi-res delta dhr
% [tmp,d] = downsample_delta(dhr,16*2); X=getPredictors(d,hrf);


model = [];
issquare = all(size(HRF) == size(HRF,1));   % if square mtx, assume convolution mtx

% make sure HRF is column vector, if single vector
if size(HRF,1) == 1, HRF = HRF'; end

if isdeltamtx(stimList)  % delta matrix
    % -------------------------------------------------------------------------------------------------
    %    * If delta function
    % -------------------------------------------------------------------------------------------------

    delta = stimList;

    if ~issquare
        for i = 1:size(delta,2)

            for j = 1:size(HRF,2)
                model(:,end+1) = conv(delta(:,i), HRF(:,j));
            end
        end
    end

else
    % -------------------------------------------------------------------------------------------------
    %    * If condition function
    % -------------------------------------------------------------------------------------------------

    for i = 1:max(stimList(:,1)) % condition function

        %delta(:,i) = (stimList == i);
        delta(:,i) = cast((stimList == i),'single'); % Changed for Matlab 7.9+ compatibility -- thanks, Liane Montana and Bruce McCandliss

        if ~issquare
            for j = 1:size(HRF,2)
                model(:,end+1) = conv(delta(:,i), HRF(:,j));
            end
        end

    end
end

% -------------------------------------------------------------------------------------------------
%    * If conv. matrix
% -------------------------------------------------------------------------------------------------

if issquare % convolution matrix
    model = HRF * delta;
end

model = model(1:size(stimList,1),:);      	% eliminate extra values

% downsample, if necessary
if ~isempty(varargin)

    dsrate = varargin{1};
    [n, k] = size(model);
    nt = n ./ dsrate;

    if nt ~= round(nt), error('Length of stimList is not evenly divisible by downsampling factor.'); end

    modeli = zeros(nt , k);

    for i = 1:size(model, 2)
        xi = 1:varargin{1}:n;       % downsample rate
        t = (1:n)';
        modeli(:, i) = interp1(t, model(:, i),xi);
    end

    model = modeli;
    %model = model(1:varargin{1}:end,:); % equivalent to resample(X,1,varargin{1},0)
end

% do not do this if you want to do nonlinear saturation!
% not necessary, i think.
%model = model - repmat(mean(model),size(model,1),1);    % center predictors

end




function isdelta = isdeltamtx(stimList)

isdelta = 0;

if all( sum(stimList == 0 | stimList == 1) == size(stimList, 1) )  % all zeros or ones
    isdelta = 1;

    if min(size(stimList)) > 1 % multiple columns
        isdelta = 1;

    end

end

end

