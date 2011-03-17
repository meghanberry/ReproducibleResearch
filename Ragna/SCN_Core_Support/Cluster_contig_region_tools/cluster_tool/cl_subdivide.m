function clout=cl_subdivide(clin,varargin);

% Usage:
% 
% [clout]=adjacent_test(clin)
% [clout]=adjacent_test(...,criterion)
% [clout]=adjacent_test(...,criterion,min_adj)
% [clout]=adjacent_test(...,criterion,min_adj,adjacency)
% % [clout]=adjacent_test(...,criterion,min_adj,adjacency,k)
% 
% Note that inputs must be specified in the correct order.
% 
% In order to ensure that all the elements of clin are in the same voxel
% space, the function passes clin to check_cl.m. If there are issues in the
% cluster, some additional command line input may be requested from the
% user.
% 
% Takes the cluster structure, clin, and subdivides clusters by requiring
% that 'core' cluster voxels have criterion (a scalar) adjacent voxels.
% Adjacent voxels are those that share at least one surface, edge, or
% corner, depending on the adjacency input.
% 
% adjacency is a string, which must be 'surface' (6 connectivity scheme),
% 'edge' (18 connectivity scheme), or 'corner' (26 connectivity scheme).
% Default is 'edge'.
% 
% criterion (must be scalar) is the number of voxels a voxel must be
% adjacent too for it to be considered a 'core' cluster voxel (all voxels
% adjacent to 'core' voxels will be included in a cluster, but those that
% are adjacent to voxels that are part of a cluster but are not themselves
% 'core' voxels will not be included in the cluster). Default is 0.
% 
% min_adj (must be scalar) is the minimum number of voxels that a voxel
% must be adjacent to for it to be included in a cluster. Voxels adjacent
% to less voxels than min_adj will recieve a value of 0 in ind. Default is
% 0.
% 
% The optional input k will require that clusters in clout have at least k
% voxels (default is 1).
% 
% 
% One of the called functions (adjacent_test.m) was refactored 10/16/06.

clin=check_cl(clin);

if isempty(varargin)
    criterion=0;
    min_adj=0;
    adjacency='edge';
    kt=1;
elseif length(varargin)<2
    criterion=varargin{1};
    min_adj=0;
    adjacency='edge';
    kt=1;
elseif length(varargin)<3
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency='edge';
    kt=1;
elseif length(varargin)<4
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency=varargin{3};
    kt=1;
else
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency=varargin{3};
    kt=varargin{4};
end

XYZ=[];
Z=[];
for k=1:length(clin)
    XYZ(:,end+1:end+size(clin(k).XYZ,2))=clin(k).XYZ;
    Z(end+1:end+length(clin(k).Z))=clin(k).Z;
end

ind=adjacent_test(XYZ,criterion,min_adj,adjacency);

k=1;
while k<=max(ind)
    if length(nonzeros(ind==k))<kt
        ind(ind==k)=0;
        ind(ind>=k)=ind(ind>k)-1;
        k=k-1;
    end
    k=k+1;
end

for k=1:max(ind)
    clout(k).XYZ=XYZ(:,ind==k);
    clout(k).XYZmm=voxelToMm(clout(k).XYZ,clin(1).M);
    clout(k).Z=Z(ind==k);
    clout(k).mm_center=center_of_mass(clout(k).XYZmm,clout(k).Z);
    clout(k).numVox=length(clout(k).Z);
    clout(k).M=clin(1).M;
    clout(k).voxSize=diag(clout(k).M(1:3,1:3));
    clout(k).name='';
    clout(k).dim=clin(1).dim;
    clout(k).threshold=clin(1).threshold;
    clout(k).title=['Cluster of ' num2str(clout(k).numVox) ' voxels, created by cl_subdivide.m with criterion=' num2str(criterion) ', min_adj=' num2str(min_adj) ', adacency=' adjacency '.'];
end