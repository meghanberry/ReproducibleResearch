V =spm_vol(DB.PP);

mask = 'Activation_thresholded.img';
vm = spm_read_vols(spm_vol(mask));
wh = find(vm);  

fprintf(1,'Loading %3.0f images: 0000',length(V));

for i = 1:length(V)
    fprintf('\b\b\b\b%04d',i);
    v = spm_read_vols(V(i));

    x(:,i) = v(wh);
end
fprintf(1,'\n')

if save_intermediate
    save parcel_images_working
end

% reduce


[pc,score,latent] = princomp(zscore(x),'econ');

% Eigenvalue plot
figure('Color','w'), set(gca,'FontSize',18),bar(latent),
xlabel('Components'),ylabel('Variance')

num_to_save = sum(latent>1);


% Component plot
nmdsfig(pc,n,nms)

% clustering

classes = clusterdata(score(:,1:3),'maxclust',10,'linkage','average');