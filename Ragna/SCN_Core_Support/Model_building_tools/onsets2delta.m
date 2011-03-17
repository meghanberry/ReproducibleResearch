function [X,delta,delta_hires,hrf] = onsets2delta(ons,TR,varargin)
% [X,delta,delta_hires,hrf] = onsets2delta(ons,TR,[Len],[hrf])
% Tor Wager, 2 / 24 / 04
%
% Builds design matrix X and delta function, given cell array of onset times in s
% One cell per condition per session, e.g., ons{1} = [24 27 29 44]';
%
% Len is optional length in s for model
%
% Limitation: can only handle two events max within the same TR
%
% Outputs:
% X is model
% delta is indicator mtx, sampled at TR
% delta_hires is indicator sampled at 16 Hz
% hrf is hemodynamic response, sampled at 16 Hz
%
% DEPRECATED -- USE ONSETS2FMRIDESIGN

disp('This function is deprecated.  use onsets2fmridesign.m');

% ----------------------------------------------
% Defaults
% ----------------------------------------------

res = 16;   % resolution, in samples per second

if length(varargin) > 0 % pre-specified length
    len = varargin{1};
    for i = 1:length(ons),
        ons{i}(round(ons{i}) > len) = [];
        
        % make column vectors if needed
        if length(ons{i}) > size(ons{i},1)
            ons{i} = ons{i}';
        end
    end
    len = len .* res;
else
    len = round(max(cat(1,ons{:})) * res);
end

cf = zeros(len,1);

if length(varargin) > 1 % custom HRF; do not norm by max
    hrf = varargin{2};
else
    hrf = spm_hrf(1/res);   % canonical hrf
    hrf = hrf ./ max(hrf);
end


% ----------------------------------------------
% BUILD DESIGN
% see also getDesign5.m
% ----------------------------------------------

for i = 1:length(ons)
    
    cf2(:,i) = cf;
    cf2(round(ons{i}*res) + 1,i) = 1;   % first TR is time 0, element 1
    cf2 = cf2(1:len,:);
end

delta_hires = cf2;


X = getPredictors(delta_hires,hrf,res*TR);
X(:,end+1) = 1;

% for HRF shape estimation
len2 = ceil(max(vertcat(ons{:})) ./ TR);
cf = zeros(len2,1);
cf2 = [];

for i = 1:length(ons)
    
    tmp = 1 + round(ons{i}./TR);
    cf2{i} = cf;
    cf2{i}(tmp) = 1;  % time 0 is first element
    
    repeats = tmp(find(diff(tmp) == 0));
    cf2{i}(repeats) = cf2{i}(repeats) + 1;
end

delta = cf2;

if length(varargin) > 0,
    for i = 1:length(ons)
        delta{i} = pad(delta{i},ceil(len./res./TR - length(delta{i})));
    end
end

return

