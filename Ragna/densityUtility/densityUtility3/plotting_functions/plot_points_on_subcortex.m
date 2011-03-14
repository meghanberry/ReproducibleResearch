function [cl,han,surfhan] = plot_points_on_subcortex(DB,name,colors,condf)
    % name can be
    % 'brainstem-thalamus', 'limbic'
    % brainstem, thalamus, hypothalamus, caudate, putamen
    % or a mask image name of your choosing
    % or already-defined clusters
    %
    % classes is
    % a vector of integers (condf) or empty
    %
    % colors is
    % a cell array with one color specification
    % or one for each integer in classes
    %
    % Examples:
    % [ind,nms,condf] = string2indicator(DB.Mode);
    % wh = [1 3];
    % condf = indic2condf(ind(:,wh)); colors = {'ro' 'b^'};
    % tor_fig; [cl,han] = plot_points_on_subcortex(DB,name,colors,condf)

    hold on;

    switch name
        case 'brainstem-thalamus'
            [cl,han{1},surfhan{1}] = show_cl(DB,'brainstem',colors,condf);
            [c,h,surfhan{end+1}] = show_cl(DB,'thalamus',colors,condf);
            cl = merge_clusters(cl,c);
            han{end+1} = h;

            [c,h,surfhan{end+1}] = show_cl(DB,'hypothalamus',colors,condf);
            cl = merge_clusters(cl,c);
            han{end+1} = h;
        case 'limbic'
            [cl,han{1},surfhan{1}] = show_cl(DB,'brainstem',colors,condf);
            [c,h,surfhan{end+1}] = show_cl(DB,'thalamus',colors,condf);
            cl = merge_clusters(cl,c);
            han{end+1} = h;

            [c,h,surfhan{end+1}] = show_cl(DB,'hypothalamus',colors,condf);
            cl = merge_clusters(cl,c);
            han{end+1} = h;

            [c,h,surfhan{end+1}] = show_cl(DB,'caudate',colors,condf);
            cl = merge_clusters(cl,c);
            han{end+1} = h;

            [c,h,surfhan{end+1}] = show_cl(DB,'amygdala',colors,condf);
            cl = merge_clusters(cl,c);
            han{end+1} = h;

            [c,h,surfhan{end+1}] = show_cl(DB,'nucleus accumbens',colors,condf);
            cl = merge_clusters(cl,c);
            han{end+1} = h;
        otherwise
            [cl,han,surfhan{1}] = show_cl(DB,name,colors,condf);

    end

    axis image; axis off; axis vis3d; view(90,10); lighting gouraud
    scn_export_papersetup(400);

    return


function [cl,han,surfhan] = show_cl(DB,name,colors,condf)
    surfhan = [];
    han = [];
    
    if isstr(name)
        switch name
            case 'brainstem'
                pname = which('spm2_brainstem.img');
                surfhan = addbrain(name);
            case 'thalamus'
                pname = which('spm2_thal.img');
                surfhan = addbrain(name);
            case 'hypothalamus'
                pname = which('spm2_hythal.img');
                surfhan = addbrain(name);
            case 'caudate'
                pname = which('spm2_caudate.img');
                surfhan = addbrain(name);
            case 'amygdala'
                load amy_clusters
                cl = amy;
                surfhan = addbrain(name);

            case 'left'
                pname = which('spm2_left.img');
                surfhan = addbrain(name);
            case 'right'
                pname = which('spm2_right.img');
                surfhan = addbrain(name);

            case 'nucleus accumbens'
                P = which('NucAccumb_clusters.mat');
                load(P)
                cl = cl(1:2);
                surfhan = addbrain(name);

            case 'putamen'
                P = which('carmack_more_clusters.mat'); load(P)
                cl = put;
                surfhan= addbrain(name);
            otherwise
                pname = name;
        end

        if ~exist('cl','var')
            cl = mask2clusters(pname);
        end
    elseif isstruct(name)
        cl = name;
    end

    % get rid of small clusters
    wh = cat(1,cl.numVox);
    wh = find(wh < 50);
    cl(wh) = [];

    han = [];
    if isempty(cl)
        disp('No clusters > 50 voxels');
    else
        cl = database2clusters(DB,cl,2);

        for i = 1:length(cl)
            h = plot_points_on_brain(cl(i).XYZmm',colors,condf(cl(i).wh_points),0);
            han = [han h];
        end

        drawnow
    end

    return

