function [p,str] = cluster_surf(varargin)
    % [suface_handle,colorchangestring] = cluster_surf(varargin)
    % Tor Wager
    % surface plot of clusters on a standard brain
    %
    % inputs, in any order:
    %   clusters structures, as created in tor_extract_rois.m
    %   cell array of colors for each cluster: {[1 0 0] [0 1 0] [0 0 1]}
    %       if number of colors specified is greater than number of clusters
    %       structures entered, n+1 and n+2 colors are overlap of 2 and overlap
    %       of all clusters, respectively.
    %   file name of mat file containing brain surface vertices and faces
    %       as created with isosurface.
    %       OR special string: 'bg' 'hipp' (hcmp,thal,amy)
    %   number of mm to plot from surface (mmdeep)
    %   optional string: 'colorscale'.  This scales colors by Z-scores of voxels
    %       if used, Z scores should be in ROW vector
    %   optional string: 'heatmap'.  This is used WITH or instead of 'colorscale',
    %   and specifies
    %   that surface colors should be heatmapped to Z-scores, rather than a single
    %   color with varying hues.
    %   optional vector: reference Z-scores range, [zmin_act zmax_act
    %                   zmax_negact zmin_negact], e.g., [0 5 -5 0], use only
    %                   with 'heatmap' option
    %                   to get refZ from clusters, try:
    % clZ = cat(2,clusters.Z); refZ = [min(clZ(clZ > 0)) max(clZ) min(clZ(clZ < 0)) min(clZ)];
    %
    % color [0 1 1] (cyan) is reserved for the overlap color btwn cluster sets.
    %
    % 'colormaps'
    %       followed by custom [colors x 3] matrices for positive colors
    %       and negative colors.
    %       matlab can create some: e.g., colormap summer, jet, etc.
    %       others can be created with colormap_tor.m
    % poscm = colormap_tor([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
    % negcm = colormap_tor([0 0 1], [0 .3 1]);  % light blue to dark blue
    % create_figure('Brain Surface'); cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, 'left');
    %
    % Special keywords for sets of surfaces are:
    % left, right, bg, limbic, cerebellum, brainstem
    %
    % Other keywords:
    % 'left'
    % 'right'
    % 'amygdala'
    % 'thalamus'
    % 'hippocampus'
    % midbrain'
    %     'caudate'
    %     'globus pallidus'
    %     'putamen'
    %     'nucleus accumbens'
    %     'hypothalamus'
    %     'cerebellum'
    %
    % see also img2surf.m
    %
    %
    % example:
    % P = 'C:\tor_scripts\3DheadUtility\canonical_brains\surf_single_subj_T1_gray.mat';
    % cluster_surf(tcl,acl,P,10,{[0 1 0] [1 0 0]},'colorscale','heatmap')
    %
    % or P = h (surface handle) to use current surface in figure, and refZ
    % cluster_surf(tcl,acl,h,[3 5 -5 -3],10,{[0 1 0] [1 0 0]},'colorscale','heatmap')
    %
    % More examples:
    % cluster_surf(cl,2,'heatmap');     % brain surface.  vertices colored @2 mm
    % cluster_surf(cl,2,'bg','heatmap');    % heatmap on basal ganglia
    % cluster_surf(cl,5,'left','heatmap');  % heatmap on left surface @5 mm
    % cluster_surf(cl,2,'right','heatmap');
    %
    % A multi-color, multi-threshold display on the cerebellum
    % colors = {[1 1 0] [1 .5 0] [1 .3 .3]};
    % tor_fig;
    % sh = cluster_surf(cl{3},colors(3),5,'cerebellum');
    % cluster_surf(cl{2},colors(2),5,sh);
    % cluster_surf(cl{1},colors(1),5,sh);
    %
    % Custom colormaps:
    % create_figure('Brain Surface'); cluster_surf(cl, 2, 'heatmap','left');


    % -------------------------------------------------------------------------
    % * set up input arguments and defaults
    % -------------------------------------------------------------------------
    mmdeep = 10;
    cscale = 0;
    heatm = 0;
    viewdeg = [135 30];

    %P = which('surf_single_subj_T1_gray.mat');
    P = which('surf_spm2_brain.mat');

    % default color maps
    poscm = colormap_tor([.7 0 0], [1 1 0]);
    negcm = colormap_tor([0 0 1], [0 .5 .5]);

    actcolors = [];
    
    clind = 1;
    for i = 1:length(varargin)

        if isempty(varargin{i})
            % ignore it
        elseif isstruct(varargin{i}), cl{clind} = varargin{i}; clind = clind+1;

        elseif iscell(varargin{i}), mycolors = varargin{i};

        elseif isstr(varargin{i})
            if strcmp(varargin{i},'colorscale'), cscale = 1;
            elseif strcmp(varargin{i},'heatmap'), heatm = 1; cscale = 1;

            elseif strcmp(varargin{i},'colormaps')
                disp('Using custom color maps.');
                poscm = varargin{i + 1}; varargin{i + 1} = [];
                negcm = varargin{i + 2};  varargin{i + 2} = [];

            elseif  strcmp(varargin{i},'left')
                P = which('surf_spm2_left.mat'); %which('surf_single_subj_grayL.mat');
                viewdeg = [90 0];
            elseif  strcmp(varargin{i},'right')
                P = which('surf_spm2_right.mat'); %which('surf_single_subj_grayR.mat');
                viewdeg = [270 0];
            elseif  strcmp(varargin{i},'hires left')
                P = which('surf_spm2_brain_left.mat'); %which('surf_single_subj_grayL.mat');
                viewdeg = [90 0];
            elseif  strcmp(varargin{i},'hires right')
                P = which('surf_spm2_brain_right.mat'); %which('surf_single_subj_grayR.mat');
                viewdeg = [270 0];

            else P = varargin{i};
            end

        elseif all(ishandle(varargin{i})) && all(varargin{i} ~= round(varargin{i}))
            disp('Found surface patch handles - plotting on existing surfaces.');
            P = varargin{i}; % handle(s) for existing surface

        elseif any(ishandle(varargin{i})) && all(varargin{i} ~= round(varargin{i}))
            disp('Found existing surface patch handles, but some are invalid.  Check handles.');
            error('Exiting')
            
        elseif length(varargin{i}) > 1,   % it's a vector
            refZ = varargin{i};
        else    % it's a number, mmdeep
            mmdeep = varargin{i};

        end

    end

    if ~exist('mycolors', 'var')
        for i = 1:length(cl)
            mycolors{i} = rand(1,3);
        end
    end

    if isempty(P)
        disp(['Cannot find: ' P]);
        P = spm_get(1,'*mat','Choose brain surface file');
    end

    if ~exist('cl', 'var')
        disp('cluster_surf.m: No clusters to plot.  Try addbrain for brain surfaces with no activation.');
        p = []; str = [];
        return
    end

    disp('cluster_surf')
    disp('___________________________________________')
    fprintf('\t%3.0f cluster structures entered\n',length(cl))
    disp('  Colors are:')
    for i = 1:length(mycolors)
        disp([' ' num2str(mycolors{i})])
    end
    if length(mycolors) > length(cl)
        disp([' overlap color is ' num2str(mycolors{length(cl)+1})])
        ovlc = ['[' num2str(mycolors{length(cl)+1}) ']'];
    else
        ovlc = '[0 1 1]';
    end

    if length(mycolors) > length(cl)+1 & length(cl) > 2
        disp([' all overlap color is ' num2str(mycolors{length(cl)+2})])
        aovlc = ['[' num2str(mycolors{length(cl)+2}) ']'];
    else
        aovlc = '[1 1 1]';
    end

    disp([' Surface stored in: ' P])

    fprintf(' Building XYZ coord list\n');

    % -------------------------------------------------------------------------
    % * build xyz list
    % -------------------------------------------------------------------------
    for i = 1:length(cl)
        xyz{i} = cat(2,cl{i}.XYZmm)';

        if cscale || heatm
            for j = 1:length(cl{i})
                if size(cl{i}(j).Z,1) > size(cl{i}(j).Z,2)
                    cl{i}(j).Z = cl{i}(j).Z';
                end
            end
            Z{i} = cat(2,cl{i}.Z)';

            % order voxels from lowest to highest, so that peak colors
            % appear because they are plotted last
            tmp = [xyz{i} Z{i}];
            tmp = sortrows(tmp,4);
            xyz{i} = tmp(:,1:3); Z{i} = tmp(:,4);

            if ~cscale % if heat map only, set mycolor{1} = [1 1 1]
                mycolors{1} = [1 1 1];
            end
        end
    end

    % -------------------------------------------------------------------------
    % * build function call
    % -------------------------------------------------------------------------
    fprintf(' Building color change function call\n');

    if length(cl) > 2 && exist('aovlc') == 1
        str = ['[c,alld] = getVertexColors(xyz{1},p,mycolors{1},[.5 .5 .5],' num2str(mmdeep) ',''ovlcolor'',' ovlc ',''allcolor'',' aovlc];
    else
        str = ['[c,alld] = getVertexColors(xyz{1},p,mycolors{1},[.5 .5 .5],' num2str(mmdeep) ',''ovlcolor'',' ovlc];
    end

    if heatm
        str = [str ',''colorscale'',actcolors{1}'];
    elseif cscale
        str = [str ',''colorscale'',Z{1}'];
    end


    for i = 2:length(cl)
        str = [str ',''vert'',xyz{' num2str(i) '},mycolors{' num2str(i) '}'];
        if heatm
            str = [str ',''colorscale'',actcolors{' num2str(i) '}'];
        elseif cscale,
            str = [str ',''colorscale'',Z{' num2str(i) '}'];
        end
    end

    str = [str ');'];


    % ------------------------------------------------------------
    % for heatmap option: get actcolors
    % -------------------------------------------------------------
    if heatm
        fprintf(' Getting heat-mapped colors\n');
        if exist('refZ') == 1
            actcolors = get_actcolors(Z, refZ, poscm, negcm);
        else
            actcolors = get_actcolors(Z, [], poscm, negcm);
        end
    end

    % -------------------------------------------------------------------------
    % * run brain surface
    % -------------------------------------------------------------------------
    if ishandle(P)      % no input file, use existing handle
        fprintf(' Using existing surface image\n');
        fprintf(' Running color change.\n');
        for i = 1:length(P)
            p = P(i);
            disp([' eval: ' str])
            eval(str)
        end
    else
        % we have either an input file or a special string ('bg')
        fprintf(' Loading surface image\n');
        [dtmp,ftmp,etmp]=fileparts(P);

        if strcmp(etmp,'.mat')

            load(P);

            %%figure
            p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
                'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
            lighting gouraud;camlight right
            axis image;
            lightRestoreSingle(gca);
            %myLight = camlight(0,0);set(myLight,'Tag','myLight');
            %set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightFollowView

            view(viewdeg(1),viewdeg(2));
            drawnow


            % -------------------------------------------------------------------------
            % * run color change
            % -------------------------------------------------------------------------
            fprintf(' Running color change.\n');
            disp([' eval: ' str])
            eval(str);


            % this for subcortex stuff
        elseif strcmp(P,'bg')
            p = [];
            myp = addbrain('caudate');p = [p myp];
            run_colorchange(myp,str,xyz,mycolors);

            myp = addbrain('globus pallidus');p = [p myp];
            run_colorchange(myp,str,xyz,mycolors);

            myp = addbrain('putamen');p = [p myp];
            run_colorchange(myp,str,xyz,mycolors);

            set(myp,'FaceAlpha',1);

            axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)


            % % %             P = which('Tal_Cau.img'); %which('ICBM_caudate.img');   %('Tal_Cau.img'); %
            % % %             figure('Color','w');[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
            % % %             if findstr(P,'Tal_Cau.img'), str(49:58) = '[.5 .6 .6]'; delete(p(1)); p = p(2);,eval(str),end
            % % %
            % % %             P = which('Tal_Glo.img');
            % % %             [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
            % % %             if findstr(P,'Tal_Glo.img'), str(49:58) = '[.5 .6 .5]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
            % % %
            % % %             P = which('Tal_Put.img');
            % % %             [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
            % % %             if findstr(P,'Tal_Put.img'), str(49:58) = '[.5 .5 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
            % % %
            % % %             [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -15],'filename','brain_render_T1');
            % % %             set(p(1),'FaceColor',[.6 .4 .3]); colormap copper;material dull;axis off
            % % %             %h = findobj('Type','Light'); delete(h);
            % % %             %[az,el]=view;lightangle(az,el); lightangle(az-180,el-60);
            % % %             [h,az,el] = lightRestoreSingle(gca);lightangle(az-180,el-60);

            % % %         elseif strcmp(P,'hipp')
            % % %             P = which('Tal_Hip.img');
            % % %             figure('Color','w');[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
            % % %             if findstr(P,'Tal_Hip.img'), str(49:58) = '[.5 .6 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
            % % %
            % % %             P = which('Tal_Amy.img');
            % % %             [p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
            % % %             if findstr(P,'Tal_Amy.img'), str(49:58) = '[.4 .4 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
            % % %
            % % %
            % % %             P = which('carmack_thal_bstem.mat');
            % % %             load(P)
            % % %             imageCluster('cluster',midbrain,'color',[.7 .3 0],'alpha',.5);
            % % %             imageCluster('cluster',thal,'color',[0 .8 .3],'alpha',.5);
            % % %
            % % %             %P = which('Tal_Tha.img');
            % % %             %[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
            % % %             %if findstr(P,'Tal_Tha.img'), str(49:58) = '[.5 .6 .5]'; eval(str); end
            % % %
            % % %             [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -20],'filename','scalped_single_subj_T1');
            % % %             set(p(1),'FaceColor',[.6 .4 .3]); colormap copper;material dull;axis off
            % % %             %h = findobj('Type','Light'); delete(h); [az,el]=view;lightangle(az,el); lightangle(az-180,el-60);
            % % %             [h,az,el] = lightRestoreSingle(gca);lightangle(az-180,el-60);
            % % %
            % % %         elseif strcmp(P,'amy')
            % % %             cl2 = load('/Users/tor/Documents/tor_scripts/3DheadUtility/Tor_Rois/amygdala_tal_mni_clusters.mat')
            % % %             amy = cl2.cl(3:4); % mni, not carmack
            % % %             imageCluster('cluster',amy,'color',[.4 .4 .6],'alpha',.8);
            % % %
            % % %             %P = which('rICBM_amygdala.img');
            % % %             %[p,outP,FV, cl, myLight] = mask2surface(P,0,[0 0 .5]);
            % % %             %if findstr(P,'rICBM_amygdala.img'), str(49:58) = '[.4 .4 .6]'; pp = p; p = pp(1); eval(str); p = pp(2);, eval(str);,end
            % % %
            % % %             [D,Ds,hdr,p,bestCoords] = tor_3d('whichcuts','z','coords',[0 0 -20],'filename','scalped_single_subj_T1');
            % % %             set(p(1),'FaceColor',[.6 .4 .3]); colormap copper;material dull;axis off
            % % %             [h,az,el] = lightRestoreSingle(gca);lightangle(az-180,el-60);
            % % %
            % % %
        elseif strcmp(P,'limbic')
            p = [];

            myp = addbrain('amygdala');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            myp = addbrain('hypothalamus');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            myp = addbrain('hippocampus');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            myp = addbrain('thalamus');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            myp = addbrain('nucleus accumbens');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);

            myp = addbrain('left');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            set(myp,'FaceAlpha',1);

            axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)

        elseif strcmp(P,'brainstem')
            p = [];

            myp = addbrain('thalamus');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            myp = addbrain('hypothalamus');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            myp = addbrain('brainstem');p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);
            %myp = addbrain('caudate');p = [p myp];
            %run_colorchange(myp,str,xyz, mycolors, actcolors);


            view(90,10); axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)

        elseif strcmp(P,'cerebellum') || strcmp(P,'amygdala') || strcmp(P,'hypothalamus') ...
                || strcmp(P,'thalamus') || strcmp(P,'midbrain') || strcmp(P,'caudate') ...
                || strcmp(P,'globus pallidus') || strcmp(P,'putamen') || strcmp(P,'nucleus accumbens') ...
                || strcmp(P,'hippocampus')
            % this uses addbrain and works with any of its keywords
            p = [];

            myp = addbrain(P);p = [p myp];
            run_colorchange(myp,str,xyz, mycolors, actcolors);

            view(90,10); axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)


        else
            error('Must input mat surf file or img file to convert to surf')
        end

    end % if ishandle

    lighting gouraud
    lightRestoreSingle(gca);
    material dull
    axis off
    set(gcf,'Color','w')
    %scn_export_papersetup(400);  % this will mess up movies!!

    disp('Finished!')
    disp('___________________________________________')

    return






function actcolor = get_actcolors(datavaluesets, refZ, poscm, negcm)

    % refZ is fixed reference range for colors; should be empty to use
    % range of data

    % ------------------------------------------------------------
    % for heatmap option: define color maps - biscale hot/cool
    % -------------------------------------------------------------

    % % % %     % color map - hot
    % % % %     % --------------------------------------------
    % % % %     %h1 = (0:1/99:1)';
    % % % %     h1 = linspace(.7,1,300)';
    % % % %     %h2 = ones(size(h1));
    % % % %     h2 = linspace(0,.8,200)';
    % % % %     h2 = [h2; linspace(.8,1,100)'];
    % % % %
    % % % %     h3 = zeros(size(h1));
    % % % %     h = [h1 h2 h3];
    % % % %     %h = [h1 h3 h3; h2 h1 h3; h2 h2 h1];
    % % % %     %h(1:75,:) = []; % take only red values to start
    % % % %     % in new matlab: h = colormap(hot(300));
    % % % %
    % % % %     % color map - winter
    % % % %     % --------------------------------------------
    % % % %     %h1 = (0:1/249:1)';
    % % % %     h1 = linspace(0,.5,250)';
    % % % %     %h2 = (1:-1/(249*2):.2)';
    % % % %     h2 = linspace(1,.7,250)';
    % % % %     h3 = zeros(size(h1));
    % % % %     hc = [h3 h1 h2];

    % log scale : transform
%     for i = 1:length(datavaluesets)
%         pos = datavaluesets{i} > 0;
%         datavaluesets{i}(pos) = log(datavaluesets{i}(pos));
%         
%         neg = datavaluesets{i} < 0;
%         datavaluesets{i}(neg) = -log(abs(datavaluesets{i}(neg)));
%     end

    actcolor = map_data_to_colormap(datavaluesets, poscm, negcm, refZ);


    % % % %     % -------------------------------------------------------------
    % % % %     % determine overall z-score range
    % % % %     % -------------------------------------------------------------
    % % % %
    % % % %     zrange = cat(2,sz{:});
    % % % %     tmp = zrange(zrange > 0);
    % % % %     tmpc = zrange(zrange < 0);
    % % % %
    % % % %     if ~isempty(tmp)
    % % % %         zrange = [min(tmp) max(tmp)];
    % % % %         if length(varargin) > 0, zrange = varargin{1}(1:2);,end % for reference range, use first 2
    % % % %         %zh = zrange(1):(zrange(2)-zrange(1))./224:zrange(2);
    % % % %
    % % % %         %zh =  linspace(zrange(1),zrange(2),299);
    % % % %         zh =  linspace(zrange(1)*zrange(1),zrange(2),299);
    % % % %         zh = round(zh*100);
    % % % %
    % % % %         if isempty(zh), zh = [1 1 0], end   % only one element?
    % % % %     end
    % % % %
    % % % %     if ~isempty(tmpc)
    % % % %         zrangec = [min(tmpc) max(tmpc)];
    % % % %         if length(varargin) > 0, zrangec = varargin{1}(3:4);,end % for reference range, use first 2
    % % % %         %zhc = zrangec(1):(zrangec(2)-zrangec(1))./249:zrangec(2);
    % % % %         zhc =  linspace(zrangec(1),zrangec(2),249);
    % % % %         zhc = round(zhc*100);
    % % % %         if isempty(zhc), zhc = [0 0 1], end   % only one element?
    % % % %     end
    % % % %
    % % % %     % -------------------------------------------------------------
    % % % %     % loop through sets of input coordinates
    % % % %     % -------------------------------------------------------------
    % % % %     % break up coords into list
    % % % %     % xyz2 = {}; indx = 1;
    % % % %     % for kk = 1:1000:size(xyz,1)
    % % % %     %     setwh{indx} = (kk:min(size(xyz,1),kk+1000))';
    % % % %     %     xyz2{indx} = xyz(setwh{indx},:);
    % % % %     %
    % % % %     %     indx = indx + 1;
    % % % %     % end
    % % % %     %
    % % % %     %
    % % % %     % for myset = 1:length(sets)
    % % % %
    % % % %
    % % % %
    % % % %
    % % % %     for i = 1:length(sz)
    % % % %
    % % % %         % -------------------------------------------------------------
    % % % %         % find color for each xyz
    % % % %         % -------------------------------------------------------------
    % % % %         clear h2,clear wh
    % % % %         myz = sz{i};
    % % % %
    % % % %         for j = 1:length(myz)
    % % % %             if myz(j) >= 0, docool = 0; else, docool = 1;, end
    % % % %
    % % % %             if docool,
    % % % %
    % % % %                 tmp = find(round(myz(j)*100) == zhc);
    % % % %                 if isempty(tmp),
    % % % %                     tmp = find((zhc-round(myz(j)*100)).^2 == min((zhc-round(myz(j)*100)).^2));
    % % % %                 end
    % % % %             else
    % % % %                 tmp = find(round(myz(j)*100) == zh);
    % % % %                 if isempty(tmp),
    % % % %                     tmp = find((zh-round(myz(j)*100)).^2 == min((zh-round(myz(j)*100)).^2));
    % % % %                 end
    % % % %             end
    % % % %
    % % % %             wh(j) = tmp(1);
    % % % %
    % % % %             if docool
    % % % %                 actcolor{i}(j,:) = hc(wh(j),:);
    % % % %             else
    % % % %                 actcolor{i}(j,:) = h(wh(j),:);
    % % % %             end
    % % % %
    % % % %         end
    % % % %
    % % % %     end
    % % % %



    % % % %     % -------------------------------------------------------------
    % % % %     % color scale bar - we must create by hand
    % % % %     % -------------------------------------------------------------
    % % % %     if length(varargin) == 0
    % % % %         try
    % % % %
    % % % %             zrange = cat(2,datavaluesets{:});
    % % % %             tmp = zrange(zrange > 0);
    % % % %             tmpc = zrange(zrange < 0);
    % % % %
    % % % %             if ~isempty(tmp)
    % % % %                 figure('Color','w'); subplot(4,1,1);hold on;
    % % % %                 zh2 = zh./100;
    % % % %
    % % % %                 for i = 2:size(h,1), fill([zh2(i-1) zh2(i-1) zh2(i) zh2(i)],[0 1 1 0],h(i,:),'EdgeColor','none');, end
    % % % %                 set(gca,'YTickLabel','');
    % % % %                 xlabel('Z-score','FontSize',14)
    % % % %
    % % % %                 docolbar = 0;
    % % % %             end
    % % % %
    % % % %             if ~isempty(tmpc)
    % % % %                 figure('Color','w'); subplot(4,1,2); hold on;
    % % % %                 zh2 = zhc./100;
    % % % %                 axis([0 .3 zh2(1) zh2(end)]),hold on
    % % % %                 for i = 1:size(h,1), plot([0 1],[zh2(i) zh2(i)],'Color',hc(i,:));, end
    % % % %                 set(gca,'XTickLabel',''); % ylabel('Z-score')
    % % % %                 h3 = get(gcf,'Position');
    % % % %                 set(gcf,'Position',[h3(1:2) h3(3)*.3 h3(4)*.5])
    % % % %                 docolbar = 0;
    % % % %             end
    % % % %
    % % % %         catch
    % % % %             figure; disp('Cannot make colorbar.  Only one voxel?')
    % % % %         end
    % % % %
    % % % %     end % if no reference Z
    return




function run_colorchange(myp,str,xyz, mycolors, actcolors)

    set(myp,'FaceAlpha',1);


    for i = 1:length(myp)
        p = myp(i);

        % get original color
        origcolor = get(p,'FaceColor');

        % color change
        eval(str);

        % change non-active back to original color
        vdat = get(p,'FaceVertexCData');
        wh = find(all(vdat == .5,2));
        vdat(wh,:) = repmat(origcolor,length(wh),1);
        set(p,'FaceVertexCData',vdat);

    end
    p = myp;
    lighting gouraud; lightRestoreSingle(gca);

    return





