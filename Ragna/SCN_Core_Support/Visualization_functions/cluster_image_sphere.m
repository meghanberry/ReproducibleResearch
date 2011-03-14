function [hpatch, cl] = cluster_image_sphere(cl, varargin)
% function [hpatch, cl] = cluster_image_sphere(cl, varargin)
%
% Images spheres at cluster centers
% Combine with addbrain and cluster_tools('connect3d') 
% or cluster_nmdsfig_glassbrain
% to get 3-D plots of connected blobs
%
% Outputs: patch handles, and new cl with sphere coordiates in XYZmm and
% XYZ
%
% Optional inputs:
% {'color', 'colors'}, mycolor = varargin{i+1};
% 'radius', myradius = varargin{i+1};
%
%
% Tor Wager, July 2007
%
% Usage:
% function [hpatch, cl] = cluster_image_sphere(cl)
% With optional arguments:
% [hpatch, cl] = cluster_image_sphere(cl, 'color', 'b', 'radius', 10)
% [hpatch, cl] = cluster_image_sphere(cl, 'color', {'r' 'g' 'b' etc}, 'radius', 10)
%
% Example: Given an MNI coordinate, plot a sphere on a brain surface
% my_mm_coord = [40, 46, 22]';
% create_figure('surface')
% cl = [];
% cl.XYZmm = my_mm_coord;
% cl.mm_center = my_mm_coord';
% V = spm_vol('brainmask.nii');
% cl.M = V.mat;
% [hpatch, cl] = cluster_image_sphere(cl, 'color', 'g', 'radius', 10)
% p = addbrain;
% set(p, 'FaceAlpha', 1);
% axis image
% view(135, 30); lighting gouraud; lightRestoreSingle; material dull;


    mycolor = 'r';
    myradius = 8;
    
    for i = 1:length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                % reserved keywords
                case {'color', 'colors'}, mycolor = varargin{i+1};

                case 'radius', myradius = varargin{i+1};

                otherwise, warning(['Unknown input string option:' varargin{i}]);
            end
        end
    end

    % get colors for each sphere
    if length(mycolor) == 1 || ~iscell(mycolor)
        mycolor = repmat({mycolor}, 1, length(cl));
    end

   newcl = convert_cl_xyz_to_sphere(cl(1), myradius);
    
    for i = 1:length(cl)

        newcl(i) = convert_cl_xyz_to_sphere(cl(i), myradius);

        hpatch(i) = image_isosurface(newcl(i), mycolor{i});


    end

    cl = newcl;
    
end





function cl = convert_cl_xyz_to_sphere(cl, myradius)

    newXYZmm = round(points_in_sphere(cl(1).mm_center, myradius)');

    cl(1).Z = ones(1, size(newXYZmm,2));
    cl(1).XYZmm = newXYZmm;
    cl(1).XYZ = mm2voxel(cl(1).XYZmm, cl(1).M, 1)';
    cl(1).numVox = size(cl(1).XYZ, 2);

end


function hpatch = image_isosurface(cl, mycolor)

    % controls padding to make sure we cover whole area
    padval = 1;
        
    % controls smoothness, etc.
    mythresh = .5;
    mysmoothbox = 3;
    mygaussstd = 1;

    mm_coords = cl(1).XYZmm;

    xyzmin = min(mm_coords') - padval;     % minus/plus for padding
    xyzmax = max(mm_coords') + padval;

    [X, Y, Z] = meshgrid(xyzmin(1):xyzmax(1), xyzmin(2):xyzmax(2), xyzmin(3):xyzmax(3));


    % construct volume data for area to image

    xvox = mm_coords(1,:) - xyzmin(1) + 1;
    yvox = mm_coords(2,:) - xyzmin(2) + 1;
    zvox = mm_coords(3,:) - xyzmin(3) + 1;

    V = zeros(size(X));

    for i = 1:size(xvox, 2)
        V(xvox(i), yvox(i), zvox(i)) = cl(1).Z(i);
    end

    % not needed if we have all mm points in sphere
    %VI = interp3(cl(1).XYZmm(1,:), cl(1).XYZmm(2,:), cl(1).XYZmm(3,:), cl(1).Z, X, Y, Z);

    V = smooth3(V, 'gaussian', mysmoothbox, mygaussstd);
    FV = isosurface(X,Y,Z,V, mythresh);

    hpatch = patch(FV);
    set(hpatch, 'EdgeColor', 'none', 'FaceColor', mycolor);

    drawnow

    %lighting gouraud

end