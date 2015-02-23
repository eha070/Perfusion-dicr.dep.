function qcc = convertFlowStagToCC(qmat) 
%function qcc = convertFlowStagToCC(qmat)  
% Converts a staggered Flow-Field to a cell-centered flow-field by
% averaging. Currently only 2D-Flow-Fields are supported
% 
% INPUT:
%   qmat - Cell-Array with staggered flow-field, numel(qmat)=3
% 
% OUTPUT:
%   qcc  - Cell-Array with cell-centered flow-field, numel(qcc)=3
%
%
%
%
%                                          (c)Constantin Heck, 23-Feb-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 
 
 
if nargin == 0 
	clc
	help(mfilename) 
	return; 
end


%prepare flow-field
qy = qmat{1};
qx = qmat{2};

%check if dimensionality is right
assert(size(qy,3)==1,'Currently only 2D-Flow-Fields are supported');

%get real m
m = size(qy) - [1,0];

%averaging from stg to cc
av =@(k) spdiags(ones(k,1)*[.5,.5],[0,1],k,k+1);
id =@(k) speye(k);

Ay = kron(id(m(2)),av(m(1)));
Ax = kron(av(m(2)),id(m(1)));


%output
qcc{1} = reshape(Ay*qy(:),m);
qcc{2} = reshape(Ax*qx(:),m);

 
 
 
end 