load meta061106; % loads datafile
DB=Meta_Setup(DB,10); % sets the kernel smooting to 10mm
DB=Meta_Select_Contrasts(DB)

%--------------------------------------------------------------------------
% Manually type this:
%--------------------------------------------------------------------------
Emotion % need to type Emotion to choose the data
1:8 % type this to pick all of it
%--------------------------------------------------------------------------

save SETUP DB
load SETUP
MC_Setup = Meta_Activation_FWE('setup', DB) % setting up the program
%--------------------------------------------------------------------------
% Manually type this:
%--------------------------------------------------------------------------
0
%--------------------------------------------------------------------------


load  MC_Info % loading the info needed for the MC simulation
Meta_Activation_FWE('mc',300) % can set the number of iteration to anything
Meta_Activation_FWE('mc',300) % can stop the simulation with Ctrl+C
                              % and run again to add more iterations
% this will do in total 600 iterations. The paper does 5000 in total. My
% results show 5300 iterations.

Meta_Activation_FWE('results',1)
spm_orthviews('Reposition',[0 0 -16])
spm_orthviews('Xhairs','off')

% plot the density map
warning off, spm_check_registration('Activation_proportion.img'), set(gcf, 'Resize', 'on'); warning on
spm_orthviews('Image','Activation_proportion.img');
spm_orthviews('Reposition',[0 0 -16]);
spm_orthviews('Xhairs','off');
colormap(hot);

% contrast indicator maps after kernel convolution
% not happy with the outcome of this - did not save an image with this
slice=38;
for i=[1 2 437]
    [V, maskdata] = iimg_read_img(DB.maskname, 1);
    xyzlist = V.xyzlist;
    maskdata = iimg_reconstruct_vols(maskdata, V);
    maskdata = meta_reconstruct_mask(wh_inmask(MC_Setup.unweighted_study_data(:,i)==1), xyzlist, V.dim(1:3));
    xyzPlotV=[];
    [x y]=find(maskdata(:,:,(slice-5:slice+5)));
    xyzPlotV = [xyzPlotV;x,y,repmat(zind,size(x))];
    xyzPlot = voxel2mm(xyzPlotV',V.mat);
    figure;
    [handles, wh_slice, my_z] = plot_points_on_slice(xyz(Contrast==i,:), 'slice', slice, 'color', [0 0 1],'close',10);
end;
