function lperfmat = localPerfusion(qmat,Fmat,h)
%function lperfmat = localPerfusion(qmat,Fmat,h)
%
% Sets up the local perfusion, which highly overestimates the perfusion
%
%
%
%
%                                      (c)Constantin Sandmann, 19-Mar-2016 
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
%------------------------------------------------------------------------- 
 
 
%voxel-volume
hd = prod(h);

%setup total inflow
Fsource = Fmat; Fsource(Fsource<0)=0;
Fbdry   = qmat{1}(1:end-1,:) + qmat{2}(:,1:end-1); 
F       =  Fbdry + Fsource;

%setup local perfusion
lperfmat = F./hd;