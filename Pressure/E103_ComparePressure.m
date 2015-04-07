% ------------------------------------------------------------------------- 
%                               E103_ComparePressure 
% 
% Compares two ways to compute the pressure:
%
% 1) qLap as the solution of 
%
%       -K/mu*\Delta p = Q
%
% 2) qDiv as the solution of 
%
%       div(-K/mu \nabla p) = Q
%
%
% Since in our case K and mu are constant, these should coincide.
% 
%                                          (c)Constantin Heck, 07-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 


clc;
clear;
close all;


%% setup Constantin's data pLap

E102_PressureByPoissonEquation

pLap = p;
qLap = q;


%% setup Erlend's data qLap
% load data
[prm,Fmat]   = settings;
basenameflow = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);
pathload = ['results/synt-createflowTPFA-' basenameflow '.mat'];
D = load(pathload);

%get pressure
pDiv = D.pmat;
qDiv = D.qmat;


%% compare the two results

figure(1);clf;
set(1,'name','Comparison of the pressure fields')

subplot(1,3,1);
imagesc(pLap);
axis image;
ti = sprintf('pLap in range [%1.2e,%1.2e]',min(pLap(:)),max(pLap(:)));
title(ti);


subplot(1,3,2);
imagesc(pLap);
axis image;
ti = sprintf('pDiv in range [%1.2e,%1.2e]',min(pDiv(:)),max(pDiv(:)));
title(ti);


subplot(1,3,3);
imagesc(pLap - pDiv);
axis image;
ti = sprintf('Difference (pLap-pDiv) in range [%1.2e,%1.2e]',min(pLap(:)-pDiv(:)),max(pLap(:)-pDiv(:)));
title(ti);
