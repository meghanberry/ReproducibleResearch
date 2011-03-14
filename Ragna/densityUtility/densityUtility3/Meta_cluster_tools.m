function varargout = Meta_cluster_tools(meth,varargin)
% varargout = Meta_cluster_tools(meth,varargin)
%
% This function contains multiple tools for working with clusters
% derived from Meta_Activation_FWE and Meta_SOM tools
%
%
% ------------------------------------------------------
% extract data and print a table for a set of meta-analysis clusters
%
% cl = make_table(cl,MC_Setup,['plot'],['successive'])
% Example:
% load Valence_Neg-Pos_Pos_clusters; load MC_Info
%
% cl1 = Meta_cluster_tools('make_table',cl{1},MC_Setup);
% cl = Meta_cluster_tools('make_table',cl,MC_Setup,'successive');
%
% ------------------------------------------------------
%
% ------------------------------------------------------
% get data for studies within rois
% studybyroi is studies activating in each cluster in cl. operator is "any" voxel in cluster counts
% [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,dat,[volInfo])
% [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,MC_Setup.unweighted_study_data,MC_Setup.volInfo)
% ------------------------------------------------------
%
% ------------------------------------------------------
% count studies by condition and plot [optional]
% [prop_by_condition,se,num_by_condition,n] = ...
% Meta_cluster_tools('count_by_condition',studybyset,MC_Setup.Xi,MC_Setup.wts,1)
%
% [prop_by_condition,se,num_by_condition,n] = ...
%   Meta_cluster_tools('count_by_condition',studybyroi,MC_Setup.Xi,MC_Setup.wts,1, ...
%   {'Right' 'Left'},MC_Setup.connames(1:5),{[1 0 0] [0 1 0] [1 0 1] [1 1 0] [0 0 1]});
%
% Xi = SOMResults.Xi(:,9:13);
% nms = SOMResults.alltasknms(9:13)
% w = SOMResults.w;
% colors = {[1 0 0] [0 1 0] [1 0 1] [1 1 0] [0 0 1]};
% [prop_by_condition,se,num_by_condition,n] = ...
%   Meta_cluster_tools('count_by_condition',studybyroi,Xi,w,1, ...
%   {'Right' 'Left'},nms,colors);
%
% ------------------------------------------------------
%
% Example:
% ------------------------------------------------------
% Run an analysis with Meta_Chisq_new, and then use these tools to get
% plots of regions. The lines below run the entire analysis.
% R = Meta_Chisq_new('compute',MC_Setup,'mask',mask);
% R = Meta_Chisq_new('write',R);
% [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,R.dat,R.volInfo);
% [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyset,R.Xi,R.w,1);
% ------------------------------------------------------

    switch meth

        case 'make_table'

            % cl = make_table(cl,MC_Setup,[doplot],[successiveflag])
            %
            doplot = 0; dosuccessive = 0;
            if any(strcmp(varargin,'successive')), dosuccessive = 1; end
            if any(strcmp(varargin,'plot')), doplot = 1; end

            varargout{1} = make_table(varargin{1},varargin{2},doplot,dosuccessive);

        case 'getdata'

            %[studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,MC_Setup.unweighted_study_data)
            %[studybyroi,studybyset] = getdata(cl,inputdata)
            if length(varargin) < 3
                [varargout{1},varargout{2},varargout{3}] = getdata(varargin{1},varargin{2});
            else
                [varargout{1},varargout{2},varargout{3}] = getdata(varargin{1},varargin{2},varargin{3});
            end

        case 'count_by_condition'
            %[prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',dat,Xi,w,doplot,[xnames],[seriesnames])
            %[prop_by_condition,se,num_by_condition,n] = count_by_condition(dat,Xi,w,doplot)

            if length(varargin) < 4, varargin{4} = 0; end
            if length(varargin) < 5, varargin{5} = []; end
            if length(varargin) < 6, varargin{6} = []; end
            if length(varargin) < 7, varargin{7} = []; end  %colors
            [varargout{1},varargout{2},varargout{3},varargout{4}] = count_by_condition(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7});

    end


    return



    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------

    % Get data within rois

    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------


function cl = make_table(cl,MC_Setup,doplot,dosuccessive)

    % uses Xi and Xinms to get differences among conditions
    % if there is no Xi field, we don't have differences among
    % conditions, so create a dummy one to get overall proportion
        if ~isfield(MC_Setup,'Xi')
            MC_Setup.Xi = ones(size(MC_Setup.wts));
            MC_Setup.Xinms = {'Act'};
        end
        
    if dosuccessive
        for i = 1:length(cl)
            disp(['Cluster cell ' num2str(i)])
            cl{i} = get_props_subfcn(cl{i},MC_Setup,doplot);
        end
    else
        cl = get_props_subfcn(cl,MC_Setup,doplot);
    end

    
    % build table function call
    if dosuccessive
        estr = 'cl = cluster_table_successive_threshold(cl,5';
    else
        estr = 'cluster_table(cl,1,0';
    end

    fnames = MC_Setup.Xinms;
    for i = 1:length(fnames), fnames{i} = [fnames{i} '_prop']; end
    nconds = length(fnames);

    for i = 1:nconds
        estr = [estr ',''' fnames{i} ''''];
    end
    estr = [estr ');'];

    % run table
    eval(estr)

    return
    
    % dependent on above:
    function cl = get_props_subfcn(cl,MC_Setup,doplot)
        disp(['getting clusters for local maxima at least 10 mm apart']);
        cl = subclusters_from_local_max(cl,10);
        cl = merge_nearby_clusters(cl,10,'recursive');
        
        % get proportion of points activating in each condition in each region
        disp('Getting studies that activated in each region.')
        [studybyroi,studybyset] = Meta_cluster_tools('getdata',cl,MC_Setup.unweighted_study_data,MC_Setup.volInfo);
        
        disp('Counting studies by condition')
        [prop_by_condition,se,num_by_condition,n] = Meta_cluster_tools('count_by_condition',studybyroi,MC_Setup.Xi,MC_Setup.wts,doplot);

        % get field names for conditions
        fnames = MC_Setup.Xinms;
        for i = 1:length(fnames), fnames{i} = [fnames{i} '_prop']; end
        nconds = length(fnames);

        fprintf(1,'Adding field to cl: %s\n',fnames{:});

        % store proportions in clusters for table printout and posterity
        for i = 1:length(cl)
            for j = 1:nconds
                cl(i).(fnames{j}) = 100 * prop_by_condition(i,j);
            end
        end

        return


function [studybyroi,studybyset, cl] = getdata(cl,inputdata,varargin)

    if length(varargin) > 0
        volInfo = varargin{1};
        %maskname = volInfo.fname;
    else
        disp('Using default mask. if your data has a different set of voxels, enter volInfo as input.')
        maskname = which('scalped_avg152T1_graymatter_smoothed.img');
        volInfo = iimg_read_img(maskname);
    end

    n_inmask_in = size(inputdata, 1);
    if n_inmask_in ~= volInfo.n_inmask
        fprintf('*****************************\nWARNING\n*****************************\n')
        fprintf('Voxels in input data set: %3.0f\nVoxels in volInfo: %3.0f\n', n_inmask_in, volInfo.n_inmask);
        fprintf('These must match!\n*****************************\n');
    end

    nrois = length(cl);
    nstudies = size(inputdata,2);

    studybyroi = false(nstudies,nrois);

    for i = 1:nrois
        [imgvec,maskvec] = iimg_clusters2indx(cl(i),volInfo);    %maskname);

        dat = inputdata(maskvec,:);
        
        cl(i).Z = sum(full(dat'));
        cl(i).Z_descrip = 'Unweighted sum of activating studies';

        studybyroi(:, i) = any(dat,1)';
        
        cl(i).activating_comparisons = studybyroi(:, i);
        
    end

    [imgvec,maskvec] = iimg_clusters2indx(cl,volInfo);   %maskname);
    dat = inputdata(maskvec,:);
    studybyset = full(any(dat)');

    studybyroi = full(studybyroi);

    return


    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------

    % Count studies in each region by condition and plot if asked for

    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------

function [prop_by_condition,se,num_by_condition,n] = count_by_condition(dat,Xi,w,doplot,varargin)

    if nargin < 4, doplot = 0; end

    [nstudies,ntasks] = size(Xi);
    if size(dat,1) ~= nstudies, dat = dat'; end
    if size(dat,1) ~= nstudies, error('data size does not match Xi'); end
    nregions = size(dat,2);

    % get stats for the entire matrix of SOMs
    [icon,ctxtxi,betas,num_by_condition,prop_by_condition] = meta_apply_contrast(dat', ...
        Xi,w,ones(1,ntasks));

    n = sum(Xi);
    n = repmat(n,nregions,1);
    se = ( (prop_by_condition .* (1-prop_by_condition) ) ./ n ).^.5;

    if doplot
        xnames = [];
        seriesnames = [];
        mycolors = [];
        if length(varargin) > 0, xnames = varargin{1}; end
        if length(varargin) > 1, seriesnames = varargin{2}; end
        if length(varargin) > 2, mycolors = varargin{3}; end

        tor_fig;
        barplot_grouped(prop_by_condition,se,xnames,seriesnames,'inputmeans','colors',mycolors);

        ylabel('Proportion of studies activating');
    end

    return








function activation_table(DB,testfield)
    % Table header

    fprintf(1,'Studies activating within %3.0f voxels of ROI\n',DB.radius);

    fprintf(1,'Study\tx\ty\tz\t%s\tN\tContrast #\tpoint index\tvox. dist\tmm dist\n',testfield);

    % for summary table
    values = DB.(testfield);
    levels = unique(values);
    w = DB.studyweight;
    w = w ./ mean(w);
    levelcnt = zeros(1,length(levels)); conwts = 0;

    for i=1:length(whcontrasts)

        whpts = find(con == whcontrasts(i));    % indices in DB.point lists for this contrast

        whcon = conindex(i);    % index in DB.contrast lists for this contrast

        xyzmm = [DB.x(whpts) DB.y(whpts) DB.z(whpts)];

        % convert to voxels for voxel distance
        % transpose so that 3-voxel lists are not reoriented (produces wrong
        % transform otherwise!)
        xyz = mm2voxel(xyzmm',VOL(1),1);

        % voxel distances
        d = distance(coord,xyz);
        whinradius = find(d <= DB.radius);  % indices relative to this pt list of points for this contrast w/i radius

        % overall DB indices of points for this contrast w/i radius
        whpts2 = whpts(whinradius);

        if isempty(whpts2),

            disp('Warning! Size or database mismatch. Contrast is in-radius, but no points within contrast are in-radius. Using closest:');,
            disp(['check: ' DB.PP(whcon,:)])
            whinradius = find(d==min(d)); whinradius = whinradius(1);
            whpts2 = whpts(whinradius);
            keyboard

        end


        % save summary information
        for j = 1:length(levels)
            levelcnt(i,j) = double(any(strcmp(levels{j},values(whpts2))));  % any values of this level for this contrast
            conwts(i,1) = w(whcon);
        end

        for j = 1:length(whpts2)
            pt = whpts2(j);
            fprintf(1,'%s\t%3.0f\t%3.0f\t%3.0f\t%s\t%3.0f\t%3.0f\t%3.0f\t%3.2f\t%3.2f\t\n', ...
                DB.study{pt},DB.x(pt),DB.y(pt),DB.z(pt),DB.(testfield){pt},N(pt),con(pt),pt,d(whinradius(j)),distance(mm',xyzmm(whinradius(j),:)));
        end


    end

    fprintf(1,'\n');




function [N,con,testfield] = check_fields(DB,testfield)
    if isfield(DB,'N'), N = DB.N;, elseif isfield(DB,'Subjects'), N = DB.Subjects; else, N = NaN*zeros(size(DB.x));, end
    if isfield(DB,'Contrast'), con = DB.Contrast;,
    else, con = NaN*zeros(size(DB.x));,
        disp('Warning!  You must have a field called DB.Contrasts for the table function to work properly.');
    end

    if ~isfield(DB,'connumbers'),
        error('No DB.connumbers field, which is required.  See Meta_Setup to create this field and set up analysis.');
    end

    % Define testfield (field to display in table)
    if isempty(testfield), try load SETUP testfield, catch, testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end, end

    if isempty(testfield), testfield = input('Cannot load testfield from SETUP.  Type name of field in DB to display: ','s'), go = 1;, end,

    while ~isfield(DB,testfield), disp(['NO field called ' testfield]);
        disp(DB), testfield = input('Type field name: ','s');
    end
    return
