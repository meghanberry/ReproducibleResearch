function mask = voxel2mask(voxels,maskdims)
% function mask = voxel2mask(voxels, x y z mask dimensions)
% voxels:
% 3 column vectors
% [i j k] = row, column, slice
% [x y z] in brain if brain is in analyze format
% (x is rows, y is columns, z is slices)
% Tor Wager, 10/17/01

    mask = (zeros(maskdims));

    for i = 1:size(voxels,1)
      mask(voxels(i,1),voxels(i,2),voxels(i,3)) = 1;
    end  

    mask = double(mask);
    
return