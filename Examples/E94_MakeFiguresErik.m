% ------------------------------------------------------------------------- 
%                               E94_MakeFigureErik 
% 
% Sets up some figures for Erik
% 
% 
% 
%                                      (c)Constantin Sandmann, 02-Mar-2016 
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
%------------------------------------------------------------------------- 

clear;
clc;
close all;

load smallDataSet.mat


%% normalization by hd

%setup some basic variables
hd = prod(prm.h);

%setup voxel-wise absolute inflow (cell-centered)
Q1 = qmat{1}(1:end-1,:);
Q2 = qmat{2}(:,1:end-1);
F  = Fmat; F(F<0)=0;

%setup true perfusion
P = (Q1+Q2+F)/hd*100*60;




%% plot results
figure(1);clf;
imagesc(P);
axis image;
colorbar('east','color','white');