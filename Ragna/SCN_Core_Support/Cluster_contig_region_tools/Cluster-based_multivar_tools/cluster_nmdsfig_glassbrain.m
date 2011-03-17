function cluster_nmdsfig_glassbrain(cl,classes,colors,sigmat,sigmat2, varargin)
    % cluster_nmdsfig_glassbrain(cl,classes,colors,sigmat,sigmat2, varargin)
    %
    % tor wager, nov. 06
    % edited: bug fix, april 2007
    %
    % classes = c.ClusterSolution.classes;
    % sigmat = c.STATS.sigmat;
    % sigmat2 = c.STATS.sigmat2;
    %
    % sigmat has signed values for significant relationships among clusters
    % colors is cell array of colors for each class (text or rgb vector)
    %
    % var args:
    % 'samefig' or 'existingfig' : do not create new figure
    % 'blobs' to image blobs
    % 'spheres' to image spheres (default)
    %
    % Example:
    % -------------------------------------------
    % Use output from mediation_direct_effects.m
    % and nmdsfig_tools.m
    %
    % sigmat = (direct_mtx > 0)' + (direct_mtx > 0);
    % cluster_nmdsfig_glassbrain(cl, c.ClusterSolution.classes, c.colors, sigmat, []);
    % i = 1; cluster_nmdsfig_glassbrain(cl, c.ClusterSolution.classes == i, c.colors(i), sigmat, []);
    % i = 2; cluster_nmdsfig_glassbrain(cl, c.ClusterSolution.classes == i, c.colors(i), sigmat, []);
    % view(135, 30); lightRestoreSingle(gca)
    % scn_export_papersetup; saveas(gcf,['Network' num2str(i)],'png');
    % 
    % Example 2:
    % ------------------------------------------
    % Output from parcel_clusters
    % create_figure('glass'); cluster_nmdsfig_glassbrain( ...
    % class_clusters{1}, ...
    % ones(length(class_clusters{1}), 1), ...
    % NMDS.basecolors(1), NMDS.stats.fdrsig, [], 'blobs');
    
    newfig = 1;
    blobstyle = 'spheres';

    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case {'samefig','existingfig'}, newfig = 0;

                case 'blobs', blobstyle = 'blobs';
                case 'spheres', blobstyle = 'spheres';


                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    fprintf('cluster_nmdsfig_glassbrain.m: ')
    fprintf('Making glass brain.\nBlob style is %s. Choices are ''blobs'' or ''spheres''\n', blobstyle);


    if newfig
        f1 = create_figure('nmdsfig_3d_glass'); 
    
    else
        f1 = findobj('Tag', 'nmdsfig_3d_glass');
        if ishandle(f1)
            figure(f1);
        else
            warning('Looking for existing figure tagged as nmdsfig_3d_glass, but can''t find it.');
            f1 = create_figure('nmdsfig_3d_glass'); 
        end
    end

    n = length(cl);

    if any(classes == 0)
        % remove these
        wh = find(classes == 0);
        cl(wh) = [];
        classes(wh) = [];

        if ~isempty(sigmat), sigmat(wh,:) = []; sigmat(:,wh) = []; end
        if ~isempty(sigmat2), sigmat2(wh,:) = []; sigmat2(:,wh) = []; end

        n = length(cl);
    end

    
    % Get colors
    for i = 1:n
        if classes(i) ~= 0
            mycolor = colors(classes(i)); mycolor = mycolor{1};
            if ischar(mycolor), mycolor = mycolor(1); end

            plotcolors{i} = mycolor;
        end
    end
        

    % Image blobs or spheres on brain
    % ----------------------------------------------------------
    switch blobstyle
        

        
        case 'blobs'

            for i = 1:n

                cltmp = cl(i);
                while cltmp.numVox < 20, cltmp = enlarge_cluster(cltmp); end

                if classes(i) ~= 0
                    h(i) = imageCluster('cluster',cltmp,'color',plotcolors{i},'alpha',1);
                end

            end

        case 'spheres'
            h = cluster_image_sphere(cl, 'colors', plotcolors);

        otherwise
            error('Unknown blobstyle argument.  blobs or spheres.');
    end



    % light lines
    % ----------------------------------------------------------
    if ~isempty(sigmat2)
        phan2 = [];
        for i = 1:n
            for j = (i+1):n
                issig = sigmat2(i,j);

                if issig && ~(sigmat(i,j))     % exclude ones that meet higher threshold

                    if sign(issig) < 0
                        out = nmdsfig_tools('connect3d',cl(i),cl(j),[.3 .3 .8], 2, .2);

                    else
                        out = nmdsfig_tools('connect3d',cl(i),cl(j),[.3 .3 .3], 2, .2);
                    end
                    phan2(end+1) = out.h;

                end
            end
        end
        drawnow
    end

    % heavy lines
    % ----------------------------------------------------------
    phan = [];
    for i = 1:n
        for j = (i+1):n
            issig = sigmat(i,j);

            if issig

                if sign(issig) < 0
                    out = nmdsfig_tools('connect3d',cl(i),cl(j), 'b', 3, .2);
                else
                    out = nmdsfig_tools('connect3d',cl(i),cl(j), 'k', 3, .2);
                end
                phan(end+1) = out.h;

            end
        end
    end
    drawnow



    % other elements
    % ----------------------------------------------------------

    bstem = addbrain('brainstem');
    cblm = addbrain('cerebellum');
    rt = addbrain('right');
    set(rt,'FaceColor',[.5 .5 .5]);
    set(rt,'FaceAlpha',.1)

    lf = addbrain('left');
    set(lf,'FaceColor',[.5 .5 .5]);
    set(lf,'FaceAlpha',.1)

    view(270,5)

    lighting gouraud
    lightRestoreSingle(gca)
    axis image; axis vis3d
    axis off
    material dull

    if input('Make movie? (1/0) ');

        mov = movie_tools('rotate',90,15,[],4);
        mov = movie_tools('rotate',0,90,mov,4);
        mov = close(mov);

    end