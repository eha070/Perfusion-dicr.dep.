function [v] = interp3c(im,c)
%INTERP3C Trilinear 3D interpolation of coordinates
%
%   V = INTERP3C(IM,C) trilinear interpolation of coordinates in the image
%   IM and for coordinates C. Returning the interpolated values
%
%   NB: C must be in voxel coordinates, assuming grid as in NDGRID
%
%   From http://local.wasp.uwa.edu.au/~pbourke/miscellaneous/interpolation/
%   Vxyz = 	V000 (1 - x) (1 - y) (1 - z) +
%   V100 x (1 - y) (1 - z) +
%   V010 (1 - x) y (1 - z) +
%   V001 (1 - x) (1 - y) z +
%   V101 x (1 - y) z +
%   V011 (1 - x) y z +
%   V110 x y (1 - z) +
%   V111 x y z 
%
%

dim = size(im);
ndim = numel(dim);
nc = size(c,1);
ndimc = size(c,2);

gpu = 0;
if isa(im,'gdouble') || isa(im,'gsingle')
    gpu = 1;
end;

if gpu
    ini = gones(nc,1);
else
    ini = ones(nc,1);
end;

% if 2D coordinates
if ndimc == 2
    c = [c ini ini];
end;
% if 3D coordinates
if ndimc == 3
    c = [c ini];
end;

% dim = size(im);
% ndim = numel(dim);
% if numel(dim) == 2
%     dim = [dim 1];
% end;

% out of bound
c(c < 1) = 1;
for i = 1 : ndim 
    c1 = c(:,i);
    ind = c1 > dim(i);
    c(ind,i) = dim(i);    
end;
clear c1;

cfl = floor(c);

chere = floor(c);
v000 = convert(im,dim,chere);

chere = [ceil(c(:,1))    floor(c(:,2))   floor(c(:,3))  floor(c(:,4))];
v100 = convert(im,dim,chere);

chere = [floor(c(:,1))   ceil(c(:,2))    floor(c(:,3))  floor(c(:,4))];
v010 = convert(im,dim,chere);

chere = [floor(c(:,1))   floor(c(:,2))   ceil(c(:,3))  floor(c(:,4))];
v001 = convert(im,dim,chere);

chere = [ceil(c(:,1))    floor(c(:,2))   ceil(c(:,3))  floor(c(:,4))];
v101 = convert(im,dim,chere);

chere = [floor(c(:,1))   ceil(c(:,2))    ceil(c(:,3))  floor(c(:,4))];
v011 = convert(im,dim,chere);

chere = [ceil(c(:,1))    ceil(c(:,2))    floor(c(:,3))  floor(c(:,4))];
v110 = convert(im,dim,chere);

chere = ceil(c);
v111 = convert(im,dim,chere);

% delta change within voxel, < 1
x = c(:,1) - cfl(:,1);
y = c(:,2) - cfl(:,2);
z = c(:,3) - cfl(:,3);

v = v000.*(1-x).*(1-y).*(1-z) + ...
    v100.*x.*(1-y).*(1-z) + ...
    v010.*(1-x).*y.*(1-z) + ...
    v001.*(1-x).*(1-y).*z + ...
    v101.*x.*(1-y).*z + ...
    v011.*(1-x).*y.*z + ...
    v110.*x.*y.*(1-z) + ...
    v111.*x.*y.*z;

%----------------------------------------

% get values
function [v] = convert(im,dim,c)

% if ndim <= 3
%     ind = sub2indfast(dim,c);
% elseif ndim == 4
ind = sub2ind(dim,c(:,1),c(:,2),c(:,3),c(:,4));
% end;
v = im(ind);

