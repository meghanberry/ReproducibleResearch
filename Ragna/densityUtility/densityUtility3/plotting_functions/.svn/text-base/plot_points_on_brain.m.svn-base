function h = plot_points_on_brain(XYZ,varargin)
    % function handles = plot_points_on_brain(XYZ,varargin)
    %
    % This function plots a 3-column vector of xyz points (coordinates from
    % studies) on a glass brain.  Four different views are created in 2
    % figures.
    %
    % An optional 2nd argument is a color, or list of colors.
    % An optional 3rd argument is a vector of integers to classify the points
    % into groups.  Each group gets a color in the colors vector.
    %
    % example:
    % plot_points_on_brain(saddecrease,{'go'});
    % plot_points_on_brain(sadpts,{'bo' 'rs' 'gv'},methi);
    % plot_points_on_brain(saddecrease,{cell vec of all colors});
    %
    % by tor wager
    % Modified aug 06 by tor wager to add option to not add brain surface
    % plot_points_on_brain(saddecrease,{'go'},[],0);
    h = [];
    addbrainsurface = 1;


    if size(XYZ, 2) ~= 3, error('XYZ must have n rows and 3 columns.'); end

    colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};

    if length(varargin) > 0, colors = varargin{1}; end
    if length(varargin) > 1, classes = varargin{2}; else classes = ones(size(XYZ,1),1); end
    if length(varargin) > 2, addbrainsurface = varargin{3};  end

    if isempty(classes), classes = ones(size(XYZ,1),1); end
    if ~iscell(colors), colors = {colors}; end
    if isempty(colors), colors = {'ro' 'go' 'bo' 'yo' 'co' 'mo' 'ko' 'r^' 'g^' 'b^' 'y^' 'c^' 'm^' 'k^'};, end

    % if we enter a vector of colors for each point, set up classes and colors
    if length(colors) == size(XYZ,1)
        classes = zeros(size(colors));
        coltmp = unique(colors);

        for i = 1:length(coltmp)
            wh = find(strcmp(colors,coltmp{i}));
            classes(wh) = i;
        end
        colors = coltmp';
    end

    if addbrainsurface
        f1 = figure('Color','w'); set(gca,'FontSize',18),hold on
    end

    for clas = 1:max(classes)

        h = plot3(XYZ(classes==clas,1),XYZ(classes==clas,2),XYZ(classes==clas,3), ...
            ['w' colors{clas}(2)],'MarkerFaceColor',colors{clas}(1),'MarkerSize',8);

    end

    if addbrainsurface

        addbrain;
        drawnow;

        % Add white filling to prevent see-through

        y = [-105 70 70 -105]; z = [-60 -60 78 78]; x = [0 0 0 0];
        hold on; hh = fill3(x,y,z,[.9 .9 .9]); set(hh,'EdgeColor','none','FaceAlpha',1)


        view(90,0)
        [az,el] = view;
        h = lightangle(az,el); set(gca,'FontSize',18)


        [hh1,hh2,hh3,hl,a1,a2,a3] = make_figure_into_orthviews;

    end

    return


