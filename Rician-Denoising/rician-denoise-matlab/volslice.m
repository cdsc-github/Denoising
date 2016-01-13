function volslice(f,slice,i)
%VOLSLICE  Plot slices from a 3D volume
%   VOLSLICE(f) plots an x slice, y slice, and z slice through the center
%   of the volume.
%
%   VOLSLICE(f,dim,i) plots the ith slice from dimension dim, where dim is
%   'x', 'y', or 'z'.
%
%   The dimension and index may be specified together, for example,
%   VOLSLICE(f,'x50') plots the 50th x slice.

% Pascal Getreuer 2009

if nargin < 3
    if nargin < 2
        clf;       
        subplot(2,2,2);
        volslice(f,'z',size(f,3)/2);
        subplot(2,2,3);
        volslice(f,'y',size(f,1)/2);  
         subplot(2,2,1);
        volslice(f,'x',size(f,2)/2);
        return;
    end
    
    i = max(1,round(str2double(slice(2:end))));
else
    i = max(1,round(i));
end

switch lower(slice(1))
    case 'x'
        image(permute(f(:,min(i,size(f,2)),:),[1,3,2])*255);
        xlabel z
        ylabel y
    case 'y'
        image(permute(f(min(i,size(f,1)),:,:),[2,3,1])*255);
        xlabel z
        ylabel x
    case 'z'
        image(f(:,:,min(i,size(f,3)))*255);
        xlabel x
        ylabel y
    otherwise
        error('Invalid input.');
end

axis image
colormap(gray(256));



