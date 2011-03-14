function STATS = xval_regression_multisubject(fit_method, Y, X, varargin)
% STATS = xval_regression_multisubject(fit_method, Y, X, varargin)
%
% CROSS-VALIDATED (JACKKNIFE) REGRESSION
% Leave-one observation out, predict outcomes for each missing obs.
%
% Y = outcome, obs x 1
% X = predictors, obs x variables
% fit_method = 'lasso' 'ols' 'ridge' or 'robust'
% optional inputs:
% 'pca', PCA data reduction first
% 'ndims', dims to save from PCA
% 'variable', retain max dims for each subject
%
% Single-level analysis: enter a single cell with Y and X data.
% rows are observations. (e.g., subjects)
%
% Multi-level analysis: enter a cell per subject with Y and X data.
% rows are observations within-subjects and should be independent for valid
% cross-validation.
%
% STATS = xval_regression_multisubject('lasso', pain_scores, data, 'pca', 'ndims', 'variable');
% SEE THE WIKI FOR EXAMPLES, ETC.
% Tor Wager, 3/17/09

pcsquash = 0;
doplssquash = 0;
ndims = 2;
cov = [];
verbose = 1;
verboseL = 0;
lassopath = '/Users/tor/Documents/matlab_code_external/machine_learning/lasso_rocha/lasso';
doplot = 1;
dochoose_ndims = 0;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case {'pca', 'pcsquash'}, pcsquash = 1;
            case {'pls', 'plssquash'}, pcsquash = 1; doplssquash = 1;
            case 'ndims', ndims = varargin{i+1};
            case {'cov', 'covs'}, cov = varargin{i+1};
            case 'lassopath', lassopath = varargin{i+1}; varargin{i + 1} = [];

            case {'noverb', 'noverbose'}, verbose = 0;
            case {'lowverb', 'loverb', 'lowverbose'}, verboseL = 1; verbose = 0;

            case 'choose_ndims', dochoose_ndims = 1; 

            case 'verbose' % default
                otherwise, warning(['Unknown input string option:' varargin{i}]);
    end
end
end

if verbose
    fprintf('xval_regression_multisubject\nVerbose mode (enter ''lowverbose'' to minimize or ''noverbose'' to turn off verbose output)\n')
end

N = length(Y); % number of subjects/datasets

if strcmp(fit_method, 'lasso'), addpath(lassopath), end

% Variable dimensions: retain max possible for each subject
choose_ndims();

if dochoose_ndims
    disp('Inner cross-validation; this could take a long time')
end

STATS.INPUTS.pcsquash = pcsquash;
STATS.INPUTS.ndims = ndims;
STATS.INPUTS.Y = Y;
%STATS.INPUTS.X = X;

include = 1:N;
clear subjfit subjbetas subjnbetas

for s = 1:N
    % ===========================================
    % DATASET (SUBJECT) LOOP
    % ===========================================

    %X{s} = dat{s}; %[pdat{s}(:, 1:10)];

    if verbose
        fprintf('Dataset %3.0f\n> -------------------------------\n ', s);
    elseif verboseL
        fprintf('%3.0f ', s);
    end

    % remove NaNs
    % ---------------------------------------------------------------------
    % X variables are done first, on the presumption that variables
    % (voxels) are many and observations are few
    nanvox{s} = any(isnan(X{s}));
    if all(nanvox{s}), error('All X variables appeared to have NaN values for one or more observations.'); end
    if verbose && sum(nanvox{s}), fprintf('Removed %3.0f X variables with NaNs\n', sum(nanvox{s})); end
    X{s}(:, nanvox{s}) = [];

    if isempty(cov)
        [wasnan{s}, Y{s}, X{s}] = nanremove(Y{s}, X{s});
    else
        [wasnan{s}, Y{s}, X{s}, cov{s}] = nanremove(Y{s}, X{s}, cov{s});
    end
    if all(wasnan{s}), error('All observations appeared to have NaN values for one or more variables.'); end
    if verbose && sum(wasnan{s}), fprintf('Removed %3.0f observations with NaNs\n ', sum(wasnan{s})); end

    Y_orig{s} = Y{s};

    clear fit
    tic

    % cross-validate: leave one out
    % ---------------------------------------------------------------------
    % initialize optional things
    wh_trace = [];
    v = [];

    fprintf('obs %3.0f', 0);

    for obs = 1:length(Y{s})

        fprintf('\b\b\b%3.0f', obs);

        % Select training/test data
        % -----------------------------
        jack_y = Y{s};
        jack_y(obs) = [];

        if ~isempty(cov)
            jdat = [X{s} cov{s}];
        else
            jdat = [X{s}];
        end

        test_jdat = jdat(obs, :);
        jdat(obs, :) = [];              % leave out the missing observation

        if dochoose_ndims
            % Inner cross-validation; this could take a long time
            % -----------------------------

            % update ndims(s), update wh_trace        
            [LASSO, wh_trace, best_ndims] = inner_xval_lasso(fit_method, jack_y, jdat, pcsquash, doplot, ndims, dochoose_ndims, doplssquash);

        end

        % Dim reduction
        % -----------------------------
        if pcsquash
            [v, jdat, test_jdat] = do_pcsquash(jack_y, jdat, test_jdat, ndims(s), doplssquash);
        end

        % Fit
        % -----------------------------
        % subjbetas{s} is a list of voxel weights for one subject, typically
        nvox = size(X{s}, 2);  % original voxels, not including covs or NaN voxels (will add in NaNs later)
        [subjbetas{s}, STATS.vox_weights(:, obs)] = do_fit(fit_method, jack_y, jdat, pcsquash, v, nvox, wh_trace);

        switch fit_method
            case {'logistic', 'logistictrain'}
                eta  = [1 test_jdat] * subjbetas{s};
                fit(obs, 1) = 1 ./ (1 + exp(-eta));
                
            otherwise
                % pred for the left-out observation
                if length(subjbetas{s}) == 1
                    % this is the 'intercept-only' model, where we have
                    % information only about the mean (cross-validated)
                    fit(obs, 1)  = 1 * subjbetas{s};
                else
                    fit(obs, 1)  = [1 test_jdat] * subjbetas{s};
                end
        end
        
    end

    toc

    % ---------------------------------------------------------------------
    % end cross-val loop
    % ---------------------------------------------------------------------

    subjfit{s} = fit;

    if pcsquash
        % add NaNs back in to preserve voxel order
        meanvoxweights = naninsert(nanvox{s}, mean(STATS.vox_weights, 2));

        STATS.mean_vox_weights(:, s) = meanvoxweights;
    end

    if verbose || verboseL

        create_figure('fitplot', 1, 2);
        plot(Y_orig{s}, 'o-');
        hold on; plot(fit, 'rx-');
        subplot(1, 2, 2);
        plot_correlation_samefig(fit, Y_orig{s});
        xlabel('Predicted value'); ylabel('Actual value')

    end

    myr = corrcoef(subjfit{s}, Y_orig{s}); rr(s) = myr(1,2);

    % These aren't standard outcome measures either; probably not very useful
    %     [corrr,t,p] = correlation('rho',subjfit{s}, Y_orig{s});
    %     rr_spearman(s) = corrr;

    %     [corrr,t,p] = correlation('irls',subjfit{s}, Y_orig{s});
    %     rr_robust(s) = corrr;

    pred_err(s, 1) = var(Y_orig{s} - subjfit{s}) .^ .5;

    % Null model leave-one-out prediction error
    % By predicting based on the mean of OTHER subjects, we've increased
    % the variance; we need to adjust
    jstat = jackknife(@mean, Y_orig{s});

    devs_from_mean_only_model = Y_orig{s} - jstat;

    devs_from_full_model = Y_orig{s} - subjfit{s};

    var_null(s) = var(devs_from_mean_only_model);
    var_full(s) = var(devs_from_full_model);

    rsq(s, 1) = 1 - var_full ./ var_null(s);

    % standard error of the improvement (or de-provement)
    %

    drawnow

end % dataset/subject loop

STATS.Y_orig = Y_orig;
if ~isempty(cov)
    STATS.note = 'covariates were added to predictive model'; %, if covs %removed in Y_orig stored here, if covs were entered';
else
    STATS.note = 'no covariates entered';
end

STATS.devs_from_mean_only_model = devs_from_mean_only_model;
STATS.devs_from_full_model = devs_from_full_model;
STATS.var_null = var_null;
STATS.var_full = var_full;

% Null model leave-one-out prediction error
% By predicting based on the mean of OTHER subjects, we've increased
% the variance
STATS.pred_err_null = var_null .^ .5;

STATS.pred_err = pred_err;
STATS.pred_err_descrip = 'Apparent error/loss: Root mean squared deviation from observed outcomes';
STATS.var_reduction = rsq;
STATS.var_reduction_descrip = 'Percent reduction in outcome variance due to model; negative values indicate added variance.';

STATS.subjfit = subjfit;
STATS.subjbetas = subjbetas;
STATS.r_each_subject = rr;
STATS.r_squared = STATS.r_each_subject .^ 2;
STATS.r_each_subject_note = 'r value: correlation between predicted and outcome values';
STATS.r_each_subject_note2 = 'may not be very interpretable, because null model r = -1.0';

if verbose
    disp(['Mean correlation is: ' num2str(mean(rr))])
    disp(['Mean proportion of original variance explained is: ' num2str(mean(rsq))])
    disp('The number above is based on reduction of original variance;');
    disp('it can be negative if predictors are not helpful because they add noise, increasing the overall variance!')
    disp(' ')
    fprintf('Null model pred. error is %3.2f, and full model is %3.2f\n', STATS.pred_err_null(1), STATS.pred_err(1));
    disp(' ')
end

% ================================
% --------------------------------
% Inline functions
%
% --------------------------------
% ================================

    function choose_ndims()
        if isstr(ndims) && strcmp(ndims, 'variable')
            if pcsquash == 0, warning('xval:ConflictingInputs', 'PC squash is off, so ndims input will not be used.'); end

            if verbose, fprintf('Variable number of dimensions: choosing: '); end

            clear ndims
            for i = 1:N
                if ~isempty(cov)
                    ndims(i) = min(size(X{i}, 1) - 2, size(X{i}, 2) + size(cov{i}, 2) - 2);
                else
                    ndims(i) = min(size(X{i}, 1) - 2, size(X{i}, 2) - 2);
                end
                if verbose, fprintf('%03d ', ndims(i)); end
            end
            if verbose, fprintf('\n'); end
        end

        if length(ndims) == 1 && N > 1, ndims = repmat(ndims, N, 1); end
    end

end % main function


% ================================
% --------------------------------
% Sub-functions
%
% --------------------------------
% ================================

function [v, jdat, test_jdat] = do_pcsquash(jack_y, jdat, test_jdat, ndims, doplssquash)
% PCA or PLS squash, returns jdat (scores) and v (weight vectors)
% and test_jdat, with applied weights
% test_jdat is used only to multiply by v to prepare for testing

if doplssquash
    %[v, jdat] = plssquash(jdat, jack_y, 'ndims', ndims, 'noplot');

    [T,P,W,Wstar,U,b,C,Bpls, v, Xhat,Yhat,R2X,R2Y] = PLS_nipals(jdat,jack_y, ndims);

    % was returning singles sometimes...
    v = double(v);
    
    % re-do jdat by taking PLS weighting, [ones jdat] * v (Bpls_star) = fit
    jdat = Yhat;
    test_jdat = [1 test_jdat] * v;

%     create_figure('tmp1'); plot(Yhat, jack_y, 'ko'); drawnow
%     xlabel('Yhat from PLS'); ylabel('Actual yhat');

else
    [v, jdat] = princomp(jdat, 'econ');
    v = v(:, 1:ndims);  % eigenvectors
    jdat = jdat(:, 1:ndims); % jdat now becomes the scores

    test_jdat = test_jdat * v; %(:, 1:ndims(s));  % get scores for missing test subj
end

end

% ================================
% ================================

function [betas, vox_weights, varargout] = do_fit(fit_method, jack_y, jdat, pcsquash, v, nvox, wh_trace)
        % subjbetas{s} is a list of voxel weights for one subject, typically
        % (in a linear model) (nvox+1) x 1, where +1 refers to the intercept parameter.
        % in the case of functional mediation, this could be nvox x tpoints.
        %
        % fits are the predicted outcome (y) values, given the model parameter estimates
        % and known information for the left-out observation.  In the functional mediation case,
        % a different method for producing fits is necessary.
out = [];

switch fit_method
    case 'ols'
        Xs = [ones(size(jdat, 1), 1) jdat];
        betas = pinv(Xs) * jack_y;

    case 'ridge'
        Xs = scale([ones(size(jdat, 1), 1) jdat]);
        betas = ridge(jack_y, Xs, 1);

        if pcsquash
            vox_weights = v * betas;
        else
            vox_weights = betas;
        end

    case 'robust'
        Xs = jdat;
        betas = robustfit(Xs, jack_y);

        if pcsquash
            vox_weights = v * betas;
        else
            vox_weights = betas;
        end

    case 'lasso'
        % Note: jdat should not have intercept; included automatically

        if size(jdat, 2) == 1, error('LASSO will not work right with only one predictor variable'); end

        out = lasso(jack_y, jdat);

        if isempty(out.beta)
            % this is the 'intercept-only' model, where we have
            % information only about the mean (cross-validated)
            % will not work with lasso as implemented though...
            wh_trace = [];
            betas = [];  %out.intercept(wh_trace);
        else
            if isempty(wh_trace), wh_trace = size(out.beta, 1); end % if empty, choose OLS 
            betas = [out.intercept(wh_trace) out.beta(wh_trace, :)]';
        end

    case 'logistic'

        betas = glmfit(jdat, [jack_y ones(size(jack_y))], 'binomial', 'link', 'logit');

        % plot apparent
% %                 eta  = [ones(size(jdat, 1), 1) jdat] * betas;
% %                 fit = 1 ./ (1 + exp(-eta));
% %                 create_figure('tmp'); 
% %                 plot(jdat, jack_y, 'ko');
% %                 vals = [-1.5:.1:1.5]';
% %                 etafit = [ones(length(vals), 1) vals] * betas;
% %                 fitcurve = 1 ./ (1 + exp(-etafit));
% %                 plot([-1.5:.1:1.5], fitcurve, 'k');
% %                 set(gca, 'XLim', [min(jdat) max(jdat)]);
% %                 etazero = -betas(1) ./ betas(2);  % point at which p(lie) in logistic model = .5
% %                 han = plot_vertical_line(etazero);
% %                 set(han, 'LineStyle', ':');
% %                 drawnow
                
% % %     case 'logisticpls'
% % %         now done in plssquash...
% % %          [T,P,W,Wstar,U,b,C,Bpls,Bpls_star,Xhat,Yhat,R2X,R2Y] = PLS_nipals(jdat,jack_y, 10);
% % %          betas was Bpls_star; [1 test_x] * Bpls_star = fit = trial score,
% % %          which is X in the logistic model below:
% % %          
% % %         Xtrain = [ones(size(jdat, 1), 1) jdat] * Bpls_star;  % this gives us the PLS-optimized features
% % %         betas = glmfit(Xtrain, [jack_y ones(size(jack_y))], 'binomial', 'link', 'logit');
% % % 
% % %         fit_eta = betas(1) + betas(2) * ([1 test_data]*Bpls_star)
% % %         fit = 1 ./ (1 + exp(-fit_eta));
        
    case 'logistictrain'
         [models] = classifierLogisticRegression( jdat, jack_y, [] );
        betas = models{3}(:, 1);  % intercept seems to be added as first predictor
        
    otherwise error('Unknown fit method')
end

if nargout > 2, varargout{1} = out; end

% Vox weights
% ------------
switch fit_method
    case {'ols', 'ridge', 'robust', 'logistic', 'logistictrain'}
        if pcsquash
            vox_weights = v(1:nvox, :) * betas(2:end);
        else
            vox_weights = betas(2:end);
        end

    case 'lasso'

        if pcsquash
            vox_weights = v(1:nvox, :) * out.beta(wh_trace, :)';
        else
            vox_weights = out.beta(wh_trace, :)';
        end

    otherwise error('Unknown fit method')
end


end % function

% ================================
% ================================

function [LASSO, wh_trace, best_ndims] = inner_xval_lasso(fit_method, jack_y, jdat, pcsquash, doplot, ndims, dochoose_ndims, doplssquash)

% Inner cross-validation; this could take a long time
% -----------------------------

% update ndims(s), update wh_trace

n_inner_obs = length(jack_y);
nvox = size(jdat, 2);  % original voxels, not including covs or NaN voxels (will add in NaNs later)

best_ndims = size(jdat, 2);  % initialize ndims value

% choose dims to test -- either variable, or single value
% --------------------------------------------------------------------
if pcsquash && dochoose_ndims
    dims_to_test = 2:ndims; % placeholder; select all features
else
    % run once, select all features
    dims_to_test = ndims;
end
         
% choose lasso shrinkage parameter with inner cross-validation loop
% -------------------------------------------------------------------
if strcmp(fit_method, 'lasso')
    if pcsquash
        % 2nd jdat(1,:)is test dat, but irrelevant here-- we just need dims
        [v, X2, tmp_jdat] = do_pcsquash(jack_y, jdat, jdat(1,:), dims_to_test(end), doplssquash);
        clear v tmp_jdat
    end
    out = lasso(jack_y, X2); % to get penalty values
    penalty_values = out.penalty;  %linspace(0, max(out.penalty), 10);  % resolution: 10

    fit = zeros(n_inner_obs, length(dims_to_test), length(penalty_values));
    
end

% Inner cross-validation
% --------------------------
for ii = 1:n_inner_obs
    
    % fprintf('\b\b\b%3.0f', obs);
    
    % Select data
    Y2 = jack_y;
    Y2(ii) = [];
    X2 = jdat;
    X2(ii, :) = [];
    test_jdat = jdat(ii, :); 
    
    % Dim reduction: with max number
    % We will use selected numbers of these
    % Must re-do PCA in inner loop to avoid bias
    % -----------------------------
    if pcsquash
        [v, X2, test_jdat] = do_pcsquash(Y2, X2, test_jdat, dims_to_test(end), doplssquash);
    end
        
    for jj = 1:length(dims_to_test)
        % Dimension selection loop
        
        % Fit
        % -----------------------------
        % we may also have a lasso shrinkage (penalty) parameter to choose
        % this gives us the whole lasso trace, though
        wh_trace = [];
        [subjbetas{ii,jj}, tmp_vox_weights, out2] = do_fit(fit_method, Y2, ...
            X2(:, 1:dims_to_test(jj)), pcsquash, v(:, 1:dims_to_test(jj)), ...
            nvox, wh_trace);

        clear tmp_vox_weights
        
        % pred for the left-out observation

        
        if strcmp(fit_method, 'lasso')
            % get trace from out struct
            % penalty values vary across validation loop iterations
            % fit = fit values for each possible choice of penalty
            % test_jdat is component scores if pcsquash, otherwise dims to
            % test is all
            fit(ii, jj, 1) = out2.intercept(end)' + test_jdat(1:dims_to_test(jj)) * out2.beta(end, :)';

            % lfit: lasso fit for each trace
            lfit = out2.intercept' + test_jdat(1:dims_to_test(jj)) * out2.beta';
            % interpolate to standard trace values
            lfit = interp1(out2.penalty, lfit, penalty_values);
            
            fit(ii, jj, 1:length(penalty_values)) = lfit;
            % this should go on 3rd dim of fit

            % trials x penalty choices
            %             fit_inner_xval{jj}(ii, :) = interp1(out2.penalty, fitiijj, penalty_values);
            %             err_inner_xval{jj}(ii, :) = jack_y(ii) - fit_inner_xval{jj}(ii, :);
        else
            
            % this is for the OLS solution (max trace value)
            if length(subjbetas{ii,jj}) == 1
                % this is the 'intercept-only' model
                fit(ii, jj, 1)  = 1 * subjbetas{ii,jj};
            else
                fit(ii, jj, 1)  = [1 test_jdat] * subjbetas{ii,jj};
            end
        end
        

        
    end % Dims loop
end % ii inner xval loop

for j = 1:length(dims_to_test)
    
    for k = 1:size(fit, 3)
        vec = squeeze(fit(:, j, k));
        cc = corrcoef(vec, jack_y);
        r(j, k) = cc(1,2);
        pe(j, k) = sqrt(sum((jack_y - vec) .^2));
    end
    
end

[wh_ndims, wh_penalty] = find(pe == min(pe(:)));
wh_ndims = wh_ndims(1);
wh_penalty = wh_penalty(1);

best_penalty = penalty_values(wh_penalty);
best_dims = dims_to_test(wh_ndims);
best_pe = pe(wh_ndims, wh_penalty);

if doplot
    create_figure('Error by dims and penalty values');

    [Xf, Yf] = meshgrid(penalty_values, dims_to_test);
    surf(Xf, Yf, pe)
    xlabel('Penalty values')
    ylabel('Dimensions')
    view(135, 30);
    hold on;
    plot3(best_penalty, best_dims, best_pe, 'ko','MarkerFaceColor', 'r', 'MarkerSize', 8);

end

LASSO = struct('best_penalty', best_penalty, 'best_dims', best_dims, 'best_pe', best_pe, ...
    'wh_ndims', wh_ndims, 'wh_penalty', wh_penalty);

end % function
