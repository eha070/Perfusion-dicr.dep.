function CBV = cbvEstimation(C,timeline,aif)
%function CBV = cbvEstimation(C,timeline,aif)
%
% Calculates CBV for a single tissue curve as 
%
%   CBV = int(C)/int(aif).
%
%
% 
% REFERENCES:
% [1] Sourbron: Classical Models in MR Perfusion
%
%
%                                          (c)Constantin Heck, 09-Aug-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


%setup main parameters
k        = numel(timeline);
timeline = timeline(:);
C        = C(:);
aif      = aif(:);

%get vector for integration
e  = timeline(2:end)-timeline(1:end-1);
Av = spdiags(1/2*ones(k,2),[0,1],k-1,k);
e  = e'*Av; e=e(:);


CBV = (e'*C)/(e'*aif(:));










end
