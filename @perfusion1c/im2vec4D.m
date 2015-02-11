function [data] = im2vec4D(varargin)

f = varargin{1};
dim = size(f);
dim3 = dim(1:3);
if nargin == 1
    roi = ones(dim3);
else
    roi = varargin{2};
end;

nvox = sum(roi(:));

% rearrange the data into vectors n x time
ntime = size(f,4);
% dim = size(f);
% dim3 = prod(dim(1:3));
data = zeros(nvox,ntime);
% save test
for i = 1 : ntime    
    datahere = f(:,:,:,i);
    datahere = datahere(roi == 1);
    data(:,i) = datahere(:);
end;
