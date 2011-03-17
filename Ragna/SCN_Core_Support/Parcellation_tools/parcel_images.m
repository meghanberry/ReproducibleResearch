function parcel_images(image_names, extract_mask, nuisance_covs)
% function parcel_images(image_names, extract_mask, nuisance_covs)
%
% parcel_images
% Principal components + anatomical parcellation of data
% Tor Wager, v. Jan 2010
%
% This function performs the following steps, in this order
% ============================================================
% - map mask to functional space
% - extract data from in-mask voxels
% - remove nuisance covariates (before assessing connectivity)
% - data reduction (pca)
% - plot cases (detect outliers)
% - separate data into a priori anatomical regions (LBPA40 hard-coded
%      right now; downloadable from web; see wiki)
%     (save label image mapped to functional space)
% - cluster voxels in each region to get parcels
% - save parcels and images
% - NMDS on the parcels to group them into "networks" (default = use rank data)
% - Visualization of the networks
% ============================================================
%
% Outputs:
% Creates and goes to a new directory: parcel_images_output
% Outputs saved to disk include
% (1) An image with unique numerical codes for each parcel 
% (2) A 'clusters' structure containing the parcels, with image data extracted and
% averaged over voxels within each parcel
%
% Inputs:
% image_names: names of images to extract data from, and to use for
% functional parcellation
% extract_mask: a mask image
% nuisance_covs: columns of a matrix to remove from the data, so that this
% subspace is not used to determine connectivity
% e.g., nuisance_covs = SPM.xX.X(:, 1:3); if these are nuisance
% covariates...

% Go to directory: output images will be written here
% ============================================================
disp('Creating directory: parcel_images_output');
mkdir('parcel_images_output');
cd('parcel_images_output');

% ============================================================
% Map mask image to image space
% ============================================================
scn_map_image(extract_mask, image_names(1,:), 'write', 'extract_mask.img');
extract_mask = fullfile(pwd, 'extract_mask.img');

% ============================================================
% Map anatomical label image to functional space and write
% ============================================================

labelimg = which('lpba40.spm5.avg152T1.label.nii');
disp(labelimg)

if ~exist(labelimg, 'file')
    cd ..
    disp('Cannot find label image: lpba40.spm5.avg152T1.label.nii')
    error('Must be on path');
end

% note: must be nearest neighbor interp to work!!
scn_map_image(labelimg, extract_mask, 'write', 'lbpa40_labels.img');
labelimg = fullfile(pwd, 'lbpa40_labels.img');  

% ============================================================
% Extract data
% ============================================================

% Extract, using single-precision array and sequential image loading to save space
% for large datasets.  'noexpand' is a bit faster and is compatible with SPM5 and above.
[dat, maskInfo] = iimg_get_data(extract_mask, image_names, 'verbose', 'single', 'noexpand');
 
% ============================================================
% Remove nuisance covs
% ============================================================
if ~isempty(nuisance_covs)
    disp('Removing nuisance covariates')
    dat = dat - nuisance_covs * pinv(nuisance_covs) * dat;
else
    disp('No nuisance covariates specified.')
end

%%
% ============================================================
% Data reduction
% ============================================================

[pc,score,latent] = princomp(zscore(dat),'econ');

% Eigenvalue plot
figure('Color','w'), set(gca,'FontSize',18),bar(latent)
xlabel('Components'),ylabel('Variance'); drawnow

num_to_save = sum(latent>1);
disp(num_to_save)
%%
% ============================================================
% Plot cases : diagnostic
% ============================================================

% Component plot of images (are there groups of points?
% this suggests clusters of images/subjects)
create_figure('nmdsfig'); nmdsfig(score(:, 1:2)); drawnow

% clustering of cases
% unusual cases (outliers) will have different class from most subjects
subj_classes = clusterdata(score(:,1:num_to_save),'maxclust',10,'linkage','average');
create_figure('nmdsfig'); nmdsfig(score(:, 1:2), 'classes', subj_classes);
title('Plot of cases: Unusual cases are possible outliers');
%%
% ============================================================
% Read all anatomical labels
% ============================================================

[labelInfo, labels] = iimg_read_img(labelimg, 2);
labels = round(labels); % small errors in resampling (?)
in_mask_labels = labels(maskInfo.wh_inmask);
all_labels = unique(in_mask_labels);

parcel_labels = zeros(maskInfo.n_inmask, 1);

%%
% ============================================================
% Cluster within anatomical regions
% ============================================================
% Loop through anatomical labels, and cluster data based on similarity
% within each region.  'region' refers to an anatomical region specified a
% priori, and 'parcel' refers to a functionally homogenous subset of voxels
% in that region.

fprintf('Parcellating %3.0f anatomical regions. 000', length(all_labels));

for i = 1:length(all_labels)
    fprintf('\b\b\b%3.0f', i);

    current_label = all_labels(i);
    wh_in_region = in_mask_labels == current_label;

    %wh_in_region_image_space = maskInfo.wh_inmask(wh_in_region);

    % view this cluster
    current_cl = iimg_indx2clusters(double(wh_in_region), maskInfo);
    cluster_orthviews(current_cl, {[1 0 0]}, 'solid');

    % select principal component weights for this region and cluster
    vox_pcs = double(pc(wh_in_region, 1:num_to_save));

    max_parcels = max(2, round(sum(wh_in_region) ./ 50));
    vox_classes = clusterdata(vox_pcs,'maxclust',max_parcels,'linkage','average');

    % assign unique value for each parcel
    % This may NOT work will all labeling schemes!!
    current_unique_id = repmat(i * 100, sum(wh_in_region), 1) + vox_classes;
    parcel_labels(wh_in_region) = current_unique_id;


    % save clusters structure for each parcel
    wh_in_region_index = find(wh_in_region);

    for j = 1:max(vox_classes)
        wh_in_parcel = wh_in_region_index(vox_classes == j);
        in_parcel = zeros(length(wh_in_region), 1);
        in_parcel(wh_in_parcel) = 1;

        tmp_cl = iimg_indx2clusters(in_parcel, maskInfo);
        % if more than 1, just choose largest for cluster structure
        % this will not match label image exactly (!)
        howmany = cat(1, tmp_cl.numVox); 
        mymax = find(howmany == max(howmany));
        mymax = mymax(1);
        
        parcel_cl{i}(j) = tmp_cl(mymax);
    end
    % visualize the parcels for this region
    cluster_orthviews(parcel_cl{i}, 'unique', 'add');
    spm_orthviews('Reposition', current_cl.mm_center)

end

% ============================================================
% Save output
% ============================================================
parcel_cl_flat = cat(2, parcel_cl{:});
%
% view all the parcels in unique colors
cluster_orthviews(parcel_cl_flat, 'unique');

% write image with unique labels
labelimg_out_name = 'parcel_labels.img';
if size(parcel_labels, 1) == 1, parcel_labels = parcel_labels'; end
iimg_reconstruct_vols(parcel_labels, maskInfo, 'outname', labelimg_out_name);
disp(['Written: ' labelimg_out_name]);

%
% extract original image data again, and average over voxels within parcel
parcel_cl_flat = tor_extract_rois(image_names, parcel_cl_flat); 

save parcel_clusters_file parcel_cl parcel_cl_flat

%%
% ============================================================
% Multivariate networks: cluster regions into networks
% ============================================================
c = [];
doranks = 1;
c = nmdsfig_tools('get_data_matrix',c, parcel_cl_flat,'timeseries',1,[],doranks);
c = nmdsfig_tools('get_correlations',c);
[c.GroupSpace,c.obs,c.implied_dissim] = shepardplot(c.D,[]);

disp('saving key info in c variable in nmds_c_structure.mat');
save nmds_c_structure c
% at least n parcels/2, and if n parcels/2 > 5, then at least 5 but up to
% nparcels/10
max_cl_to_test = min(max(5, round(length(parcel_cl_flat) ./ 10)), round(length(parcel_cl_flat)./2));

c = nmdsfig_tools('cluster_solution',c, c.GroupSpace, 2:max_cl_to_test, 1000, []);


c.colors = cluster_orthviews_classes(parcel_cl_flat,c.ClusterSolution.classes, [], [], 1);
saveas(gcf,'network_montage','fig');
scn_export_papersetup(550);
saveas(gcf,'network_montage','png');

% apply and make network plot
c = nmdsfig_tools('apply_clusters',c);

saveas(gcf,'network_nmds_fig','fig');
scn_export_papersetup(450);
saveas(gcf,'network_nmds_fig','png');

%nmdsfig_tools('nmdsfig_plot',c, 0, 0, 'fill');

disp('saving updated key info in nmds_c_structure.mat');
save nmds_c_structure c


% name networks
c.APPLY_CLUSTER.names = cell(1, length(c.APPLY_CLUSTER.classes));
for i = 1:length(c.APPLY_CLUSTER.classes)
    wh = find(c.ClusterSolution.classes == i);
    cluster_orthviews(parcel_cl_flat(wh), c.colors(i), 'solid');
    c.APPLY_CLUSTER.names{i} = input(sprintf('Name network %3.0f: ', i), 's');
end

disp('saving updated key info in nmds_c_structure.mat');
save nmds_c_structure c

% re-make to leave open for display
cluster_orthviews_classes(parcel_cl_flat,c.ClusterSolution.classes, [], [], 1);

cd('..');
%%
end  % end function
