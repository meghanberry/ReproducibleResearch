function [ind]=adjacent_test(XYZ,varargin)

% USAGE:
%
% [ind]=adjacent_test(XYZ)
% [ind]=adjacent_test(...,criterion)
% [ind]=adjacent_test(...,criterion,min_adj)
% [ind]=adjacent_test(...,criterion,min_adj,adjacency)
%
% Note that inputs must be specified in the correct order.
%
% Returns a vector of same length as the non-tripleton dimension of XYZ,
% with numeric values indicating which cluster each voxel belongs too.
% Thus, assuming each column refers to a voxel coordinate,
% XYZ(:,ind==1) is all of the voxels in one cluster, and
% XYZ(:,ind==2) is another.
%
% adjacency is a string, which must be 'surface' (6 connectivity scheme),
% 'edge' (18 connectivity scheme), or 'corner' (26 connectivity scheme).
% Default is 'edge'.
%
% criterion (must be scalar) is the number of voxels a voxel must be
% adjacent too for it to be considered a 'core' cluster voxel (all voxels
% adjacent to 'core' voxels will be included in a cluster, but those that
% are adjacent to voxels that are part of a cluster but are not themselves
% 'core' voxels will not be included in the cluster). Setting criteria
% above the connectivity value associated with a given adjacency criteria
% will result in no clusters. Default is 0.
%
% min_adj (must be scalar) is the minimum number of voxels that a voxel
% must be adjacent to for it to be included in a cluster. Voxels adjacent
% to less voxels than min_adj will recieve a value of 0 in ind. Setting
% criteria above the connectivity value associated with a given adjacency
% criteria will result in no clusters. Default is 0.
%
%
% Note that:
% ind=adjacent_test(XYZ,'edge'); %or ind=adjacent_test(XYZ,'edge',0,0);
% ind2=spm_clusters(XYZ);
% isequal(ind,ind2)
%
% ans =
%
%      1
%
% 
% 
% Refactored 10/16/06. Script now a couple orders of magnitude faster, due
% to improvements in the loop for checking which voxels are adjacent to
% each voxel. Large (>10000 voxels) clusters should no longer take hours to
% subdivide. Also added some verbosity, along with progress reporting for
% every 1000 voxels in clusters of at least that size.

if isempty(varargin)
    criterion=0;
    min_adj=0;
    adjacency='edge';
elseif length(varargin)<2
    criterion=varargin{1};
    min_adj=0;
    adjacency='edge';
elseif length(varargin)<3
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency='edge';
else
    criterion=varargin{1};
    min_adj=varargin{2};
    adjacency=varargin{3};
end

if ~strcmp(adjacency,'corner')&&~strcmp(adjacency,'edge')&&~strcmp(adjacency,'surface')
    error('Adjacency input incorrectly specified')
end

if size(XYZ,1)~=3
    if size(XYZ,2)~=3
        error('XYZ matrix must have a tripleton dimension!')
    else
        XYZ=XYZ';
    end
end

if ~isequal(XYZ,round(XYZ))||any(XYZ(:)<1)||size(XYZ,2)~=length(unique(XYZ','rows'))
    error('XYZ input contains invalid or duplicate coordinates! If you have an affine matrix (M), use:  XYZ=mmToVoxel(voxelToMm(XYZ,M),M,''valid'')')
end

if strcmp(adjacency,'edge')&&~criterion&&~min_adj
    try ind=spm_clusters(XYZ);return,catch end
end

ind=zeros(1,size(XYZ,2));

disp('Checking for adjacent clusters')

for k=1:length(ind)
    
    if ~rem(k,1000)
        disp([int2str((k/length(ind))*100) '% complete'])
    end
    
    adjacent=zeros(3,length(ind));
    for dim=1:3
        a=((XYZ(dim,k)==XYZ(dim,:)+1)|(XYZ(dim,k)==XYZ(dim,:)-1));
        b=(XYZ(dim,k)==XYZ(dim,:));
        if ~isempty(a)
            adjacent(dim,a)=1;
        end
        if ~isempty(b)
            adjacent(dim,b)=2;
        end
    end

    [r,c]=find(adjacent);
    [r,d]=find(~adjacent(:,c));
    c(d)=[];
    c=unique(c);
    
    adj=[];
    for m=1:length(c)
        if strcmp(adjacency,'corner')
            adj=[adj c(m)];
        else
            matched=(adjacent(:,c(m))==2);
            if strcmp(adjacency,'edge')&&any(matched)
                adj=[adj c(m)];
            elseif length(nonzeros(matched))>1
                adj=[adj c(m)];
            end
        end
    end

    if length(adj)<min_adj+1||length(adj)<criterion+1
        continue
    end

    if any(ind(adj))
        ind(adj)=min(nonzeros(ind(adj)));
    else
        ind(k)=max(ind)+1;
    end
end

disp('Done checking for adjacent clusters')

k=1;
while k<=max(ind)
    if ~any(ind==k)
        ind(ind>k)=ind(ind>k)-1;
        k=k-1;
    end
    k=k+1;
end

disp('Attaching stray voxels to existing clusters')

lost=[];
for k=1:length(ind)
    
    if ~rem(k,1000)
        disp([int2str((k/length(ind))*100) '% complete'])
    end
    
    if ~ind(k)
        adjacent=zeros(3,length(ind));
        for dim=1:3
            a=((XYZ(dim,k)==XYZ(dim,:)+1)|(XYZ(dim,k)==XYZ(dim,:)-1));
            b=(XYZ(dim,k)==XYZ(dim,:));
            if ~isempty(a)
                adjacent(dim,a)=1;
            end
            if ~isempty(b)
                adjacent(dim,b)=2;
            end
        end
        
        [r,c]=find(adjacent);
        [r,d]=find(~adjacent(:,c));
        c(d)=[];
        c=unique(c);
        
        adj=[];
        for m=1:length(c)
            if strcmp(adjacency,'corner')
                adj=[adj c(m)];
            else
                matched=(adjacent(:,c(m))==2);
                if strcmp(adjacency,'edge')&&any(matched)
                    adj=[adj c(m)];
                elseif length(nonzeros(matched))>1
                    adj=[adj c(m)];
                end
            end
        end

        if length(adj)<min_adj+1
            continue
        end

        if any(ind(adj))
            if min(nonzeros(ind(adj)))==max(ind(adj))
                ind(k)=max(ind(adj));
            else
                lost(end+1)=k;
            end
        end
    end
end

disp('done')

for k=1:length(lost)
    clear num l
    adjacent=zeros(3,length(ind));
    for dim=1:3
        a=find(XYZ(dim,lost(k))==XYZ(dim,:)+1);
        a=[a find(XYZ(dim,lost(k))==XYZ(dim,:)-1)];
        b=find(XYZ(dim,lost(k))==XYZ(dim,:));
        if ~isempty(a)
            adjacent(dim,a)=1;
        end
        if ~isempty(b)
            adjacent(dim,b)=2;
        end
    end

    [r,c]=find(adjacent);
    adj=[];
    for m=1:max(c)
        if length(find(m==c))==3
            if strcmp(adjacency,'corner')
                adj=[adj m];
            else
                matched=find(adjacent(:,c(m==c))==2);
                if strcmp(adjacency,'edge')&&~isempty(matched)
                    adj=[adj m];
                elseif length(matched)>1
                    adj=[adj m];
                end
            end
        end
    end

    poss=min(nonzeros(ind(adj))):max(ind(adj));
    for m=1:length(poss)
        num(m)=length(find(ind(adj))==poss(m));
    end
    if length(find(num==max(num)))==1
        ind(lost(k))=poss(num==max(num));
    else
        a=find(num==max(num));
        for m=length(a)
            l(m)=length(find(ind==poss(a(m))));
        end
        if length(find(l==max(l)))==1
            warning(['Voxel at ' num2str(XYZ(1,lost(k))) ',' num2str(XYZ(2,lost(k))) ',' num2str(XYZ(3,lost(k))) ' shares an equal number of adjacent voxels with more than one seperate cluster. It is being associated with the largest of them.']);
            ind(lost(k))=poss(a(l==max(l)));
        else
            warning(['Voxel at ' num2str(XYZ(1,lost(k))) ',' num2str(XYZ(2,lost(k))) ',' num2str(XYZ(3,lost(k))) ' shares an equal number of adjacent voxels with more than one seperate cluster of identical size. It is being randomly associated with one of them.']);
            r=randperm(length(a));
            ind(lost(k))=poss(a(r(1)));
        end
    end
end



