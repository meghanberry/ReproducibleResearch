function [within_ste,stecons] = barplot_get_within_ste(dat,cons,covs)
% [within_ste,stecons] = barplot_get_within_ste(dat,cons,covs)
%
% Compute within-subjects standard errors based on contrasts of interest
% Adjust for between-subjects covariates
%
% Useful for barplots of conditions
%
% Example:
% data is 30 x 3
% cons = [-1 1 0; 0 -1 1]; % differences of interest within-ss
% covs = 30 x 2 covariates
% [within_ste,stecons] = barplot_get_within_ste(dat,cons,covs)

if isempty(covs), covs = ones(size(dat, 1), 1); end
    
%% Check sizes and stuff
 n = size(covs,1); 
 if n ~= size(dat,1), error('Covs and dat must have same no. of observations.'); end
 
 k = size(cons,1);
 if k ~= size(dat,2), error('Contrasts must have same no. of columns as data.'); end
 
 
%% Prepare between-subjects covariates
% remove constant, if there is one
wh = (~any(covs - repmat(mean(covs),n,1)));
covs(:,wh) = [];

% center
covs = covs - repmat(mean(covs),n,1);


%% Apply within-subjects contrasts
 
 datw = dat * cons;
 
%% Remove effects of covs
 
% add intercept
    X = [ones(n,1) (covs - repmat(mean(covs), n, 1))];
    betas = pinv(X) * datw;
    r = datw - X * betas;
    
%% within-subjects sterr, controlling for covs
    
stecons = std(r) ./ sqrt(n);
within_ste = mean(stecons);

return
