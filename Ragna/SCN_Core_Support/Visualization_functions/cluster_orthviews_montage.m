% slices_fig_h = cluster_orthviews_montage(spacing, myview, [overlay], [xhairs])
%
% cluster_orthviews_montage(cl, 'coronal');
% cluster_orthviews_montage(cl, 'sagittal');
% cluster_orthviews_montage(cl, 'axial');
%
% tor wager, aug 2006
% used in cluster_orthviews_classes

function slices_fig_h = cluster_orthviews_showcenters(spacing, myview, varargin)
    
    overlay = [];
    xhairs = 0;
    if length(varargin) > 0
        overlay = varargin{1}; end
    
    if isempty(overlay), overlay = which('scalped_single_subj_T1.img'); end
    
    whichorth = [];
    for i = 1:length(varargin)
       if ischar(varargin{i})
            switch varargin{i}

                case 'whichorth', whichorth = varargin{i+1};
                
                otherwise, warning('scn_tools:badInput', ['Unknown input string option:' varargin{i}]);
            end
        end
    end 
    
    % set up orthviews to be OK
    try % will crash if no blobs
        scn_export_spm_window('setup', overlay);
    catch
    end
    
%     set(get(findobj('Tag','Graphics'), 'Children'), 'Color', 'white');
    set(findobj('Type', 'axes', 'Parent', findobj('Tag','Graphics')), 'Color', 'white');
    if xhairs, spm_orthviews('Xhairs', 'on'); end

    % Set window to be equal/constant
    [volInf_tmp, dat] = iimg_read_img(overlay, 2);
    dat = dat(volInf_tmp.wh_inmask);
    spm_orthviews('Window', whichorth, [prctile(dat, 2) prctile(dat, 98)]);

    myviews = {'axial' 'coronal' 'sagittal'};   % for selecting SPM window
    whview = find(strcmp(myviews, myview));
    if isempty(whview), error('myview must be axial, coronal, or sagittal.'); end
    
    switch whview
        case 1
            cen = [-45:spacing:70]';
            cen = [zeros(length(cen), 2) cen];
        case 2
            cen = [-100:spacing:65]';
            cen = [zeros(length(cen),1) cen zeros(length(cen),1)];
        case 3
            cen = [-70:spacing:70]';
            cen = [ cen zeros(length(cen), 2)];
    end

   
    
    myviews2 = {'sagittal' 'coronal' 'axial' };  % for selecting coord
    whcoord = strmatch(myview, myviews2) ;

    % get text string base
    mystr = {'x = ' 'y = ' 'z = '};
    textbase = mystr{whcoord};

    % get optimal number of axes
    num_axes = size(cen, 1);
    rc = ceil(sqrt(num_axes));

    if isempty(whichorth)
        axh = get_orth_axishandles;
    else
        axh = get_orth_axishandles_whichorth(whichorth);
    end

    axh = axh(whview);

    slices_fig_h = figure; %create_figure(myview);
    set(slices_fig_h, 'Color', 'k');
    
    for i = 1:num_axes
        newax(i) = subplot(rc, rc, i);
        axis off;
     end


    for i = 1:num_axes
        spm_orthviews('Reposition', cen(i, :));

        copyobj(get(axh, 'Children'), newax(i));
        axes(newax(i));
        axis image

        h = findobj(newax(i), 'Type', 'text');
        if length(h) > 1, h= h(1); end % kludgy fix for multiple matches ***
        
            % try to set a reasonable font size
            pos = get(newax(i), 'Position');
            height = pos(3);
            fs = round(height * 70);
            set(h, 'FontSize', fs)
            % set position of text (move down and right)
            pos = get(h, 'Position');
            pos(2) = 0; %pos(2) - fs./2;
            pos(1) = 0;
            set(h, 'Position', pos)
            % set text string based on cen
            set(h, 'String', [textbase num2str(cen(i))]);
%         elseif ishandle(h)
%             delete(h);
%         end
    end

    % try to set a reasonable enlargement factor
    n = 1 + .15 * log(num_axes ./ 2);
    n = max(n, 1); n = min(n, 3);
    enlarge_axes(gcf, n)

    % set background color to print as whatever it is on screen
    set(gcf,'InvertHardcopy', 'off');
end




function axish = get_orth_axishandles

    % Get figure and axis handles
    fh = findobj('Tag', 'Graphics');
    ch = get(fh, 'Children');
    for i= 1:length(ch)
        mytype = get(ch(i), 'Type');
        wh(i) = strcmp(mytype, 'axes');
    end
    axish = ch(find(wh));

    if isempty(axish), error('SPM figure orthviews do not exist'); end

    % get which axis is which
    for i = 1:length(axish)
        poss(i, :) = get(axish(i), 'Position');
    end

    % get rid of extra axes we may have created in the 4th quadrant
    other_axes = find(any(poss(:, 1) > .45 & poss(:, 2) < .2, 2));
    axish(other_axes) = [];
    poss(other_axes, :) = [];

    % sort into order:  axial, coronal, saggital
    ssum = sum(poss(:, 1:2), 2);
    [ssum, ind] = sort(ssum);
    axish = axish(ind);

end

function axish = get_orth_axishandles_whichorth(whichorth)
% Get the axis handles for the current orthviews
global st
for i = 1:length(st.vols), wh(i) = ~isempty(st.vols{i}); end
wh = find(wh); wh = wh(whichorth);
axish = cat(1, st.vols{wh}.ax{:});
axish = sort(cat(1, axish(:).ax));

end

