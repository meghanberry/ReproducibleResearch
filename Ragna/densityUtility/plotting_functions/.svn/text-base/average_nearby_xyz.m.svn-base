function [xyzall,codesall,orderall] = average_nearby_xyz(XYZ,thresh,codes)
% [xyzall,codesall,orderall] = average_nearby_xyz(XYZ,thresh,codes)
%
% for plotting peak coordinates on the brain, for example, eliminating
% redundant nearby coordinates from the same study

u = unique(codes);


codesall = [];
xyzall = [];
orderall = [];

for i = 1:length(u)
    
    studywh = find(strcmp(codes,u{i}));
    wh = XYZ(studywh,:);
    
    dis = squareform(pdist(wh));

    dis = dis < thresh;
    
    %howmany = sum(dis,2);
    %[dummy,order] = sort(1./howmany);
    %howmany = howmany(order,:);
    %dis = dis(order,:);
    
    [u2,order2] = unique(dis,'rows');
    
    xyzout = [];
    for j = 1:size(u2,1)
        
        tmpxyz = wh(find(u2(j,:) == 1),:);
        xyzout(j,:) = mean(tmpxyz,1);
        
    end
    
    xyzall = [xyzall; xyzout];
    codesall = [codesall; repmat(u(i),size(xyzout,1),1)];
    orderall = [orderall; studywh(order2)];
end



return