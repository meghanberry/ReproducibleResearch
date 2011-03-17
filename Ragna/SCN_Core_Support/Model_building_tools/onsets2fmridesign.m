% [X, delta, delta_hires, hrf] = onsets2fmridesign(onsets, TR, [len], [custom hrf or basis set name],[other optional args, i.e., norm])
% Tor Wager, 2 / 24 / 04
% Modified/updated 2/08 by Tor to add capacity for durations
%
% Builds design matrix X and delta function, given cell array of onset times in s
% One cell per condition per session, e.g., ons{1} = [24 27 29 44]';
%
% Inputs:
%   onsets:
%   - 1st column is onsets for events,
%   - 2nd column is optional durations for each event
%   - Enter single condition or cell vector with cells for each condition (each event type).
%   TR: RT in seconds
%   len: optional length in s for model, or [] to skip if using additional
%   args
%   hrf name: a string used by spm_get_bf.m, or [] to skip
%   'norm' : mean-center, orthogonalize, and L2-norm basis set
%
% Limitation: can only handle two events max within the same TR
%
% Outputs:
%   X: model
%   delta: indicator mtx, sampled at TR
%   delta_hires: indicator sampled at 16 Hz
%   hrf: hemodynamic response function, sampled at 16 Hz
%
% X = onsets2fmridesign(ons, TR, size(imgs, 1) .* TR, 'hrf (with time derivative)');


function [X, delta, delta_hires, hrf] = onsets2fmridesign(ons, TR, varargin)
    % ----------------------------------------------
    % Defaults
    % ----------------------------------------------

    res = 16;   % resolution, in samples per second

    % Other optional inputs after fixed ones
    % - - - - - - - - - - - - - - - - - - - -
    donorm = 0;

    for i = 3:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}

                % reserved keywords
                case 'norm', donorm = 1;
               
                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end
    
    
    if ~iscell(ons), ons = {ons}; end

    % Optional: pre-specified length, or []
    % - - - - - - - - - - - - - - - - - - - -
    
    if ~isempty(varargin) && ~isempty(varargin{1}) % pre-specified length
        
        len = ceil(varargin{1});
        for i = 1:length(ons)
            % make column vectors if needed
            if length(ons{i}) > size(ons{i}, 1)
                disp('Warning! onsets2fmridesign thinks you''ve entered row vectors for onsets and is transposing.  Enter col vectors!');
                ons{i} = ons{i}';
            end

            % make sure nothing goes past end of session
            past_end = round(ons{i}(:,1)) > len;
            if any(past_end)
                fprintf('Warning! onsets2fmridesign found %d onsets after the end of the session and is omitting them.\n', sum(past_end)); 
            end
            ons{i}(past_end,:) = [];
        end
        len = ceil(len .* res);
    else
        len = round(max(cat(1, ons{:})) * res);
        len = ceil(len/(res*TR)) * res*TR;
        len = len(1);
    end

    cf = zeros(len, 1);

     % Optional: basis set name, or []
    % - - - - - - - - - - - - - - - - - - - -   
    if length(varargin) > 1 
        if ischar(varargin{2})
        % basis set: Any of the named basis sets from SPM  spm_get_bf

            % bf = spm_get_bf(struct('name', 'hrf (with time and dispersion derivatives)', 'length', 30, 'dt', 1));
            bf = spm_get_bf(struct('name', varargin{2}, 'length', 30, 'dt', 1/res));
            hrf = bf.bf;

        else
            % custom HRF; do not norm by max
            hrf = varargin{2};
        end
    else
        hrf = spm_hrf(1/res);   % canonical hrf
        hrf = hrf ./ max(hrf);
    end

    % ----------------------------------------------
    % Normalize basis set, if requested
    % ----------------------------------------------
    
    if donorm
        % Mean-center
        hrf = scale(hrf, 1);
        kk = size(hrf, 2);

        % Orthogonalize "later" basis funtions wrt "earlier" (assumed to be more important) ones
        if size(hrf, 2) > 1, hrf = spm_orth(hrf); end

        % Normalize to L2 norm
        for ii = 1:kk
            hrf(:, ii) = hrf(:, ii) ./ norm(hrf(:, ii));
        end
    end
    
    % ----------------------------------------------
    % BUILD DESIGN
    % see also getDesign5.m
    % ----------------------------------------------

    for i = 1:length(ons)

        cf2(:,i) = cf;                      % add empty column for this condition
        cf2(round(ons{i}(:,1)*res) + 1, i) = 1;   % specify indicators; first TR is time 0, element 1

        we_have_durations = size(ons{i}, 2) > 1;
        if we_have_durations
            for j = 1:size(ons{i}, 1)
                dur_in_s = ons{i}(j, 2);
                first_sample = round(ons{i}(j, 1)*res) + 1;
                end_in_samples = min(len, round(first_sample + dur_in_s * res));

                cf2(first_sample:end_in_samples, i) = 1; % add to condition function
            end
        end

        cf2 = cf2(1:len,:);
    end

    delta_hires = cf2;


    X = getPredictors(delta_hires, hrf, res*TR);
    X(:,end+1) = 1;

    % for HRF shape estimation
    len2 = ceil(max(vertcat(ons{:})) ./ TR);
    cf = zeros(len2(1), 1);
    cf2 = [];

    for i = 1:length(ons)
        tmp = 1 + round(ons{i}(:,1)./TR);
        cf2{i} = cf;
        cf2{i}(tmp) = 1;  % time 0 is first element

        repeats = tmp(find(diff(tmp) == 0));
        cf2{i}(repeats) = cf2{i}(repeats) + 1;
    end

    delta = cf2;

    if ~isempty(varargin)
        for i = 1:length(ons)
            delta{i} = pad(delta{i}, ceil(len./res./TR - length(delta{i})));
        end
    end
end

