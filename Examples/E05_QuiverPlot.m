% ------------------------------------------------------------------------- 
%                               E05_QuiverPlot 
% 
% Loads qmat and shows a quiver plot of the flow
% 
% 
% 
%                                          (c)Constantin Heck, 23-Feb-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 




clear;
clc;
close all;

load('synt-createflowTPFA-phi-flat-K-flat-dim-128.mat')


%%
%setup domain
omega = [0,1,0,1];
m     = size(pmat);


%% convert flow

%get flow
qcc = perfusion1c.convertFlowStagToCC(qmat);
qy  = qcc{1};
qx  = qcc{2};


%get meshgrid to plot qcc, the cell-centered flow
ygv   = linspace(omega(2),omega(1),m(1)); %from (0 to 1)
xgv   = linspace(omega(3),omega(4),m(2)); %this goes from (1 to 0) (matlab, duh)
[x,y] = meshgrid(xgv,ygv);


%create quiver plot. Notice that flow in x-direction is reversed
quiver(y(:),x(:),qy(:),-qx(:),2)