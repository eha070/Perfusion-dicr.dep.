% ------------------------------------------------------------------------- 
%                               E98_TestPressure
% 
% Checks, if p fulfills the relationship -K/mu \Delta p = Q.
% To do this, the laplacian is calculated
% 
% 
% 
%                                          (c)Constantin Heck, 06-Apr-2015 
%                                                http://mic.uni-luebeck.de
% ------------------------------------------------------------------------- 

clear
close all;
clc


%% prepare data
% load data
[prm,Fmat]   = settings;
basenameflow = perfusion1c.providenameflow(prm.phiopt,prm.Kopt,prm.dim);
pathload = ['results/synt-createflowTPFA-' basenameflow '.mat'];
D = load(pathload);
pathload = './results/synt-createindicatorpde-phi-flat-K-flat-dim-64-aif-gamma-T-90-red-25.mat';
E = load(pathload);

%get physiological dimensions
m        = prm.dim(1:2);
h        = prm.h;
voxelvol = prod(prm.h);
dim      = numel(h);

%get Cmat, qmat and phimat
C    = E.Cmat;
AIF  = reshape(E.aifval,[],1);
qmat = D.qmat;
phi  = D.phimat;
k    = numel(AIF);

%% do the testing: 
% We know that dC/dt = div(q*c) + FIn*AIF - FOut*c
% reformulation: dC/dt = div(q)*c + q*GRADc + FIn*AIF - FOut*c
% We will calculate:
% dCdt = divqc + qGRADc + FInAIF - FOutc
%
%REMARK: This is probably not the best idea: Maybe I should just sum up
%differences over the faces


%setup gradient operator with Neumann-bc
e  =@(i) ones(m(i),1);
id =@(i) speye(m(i));

D1 = spdiags(e(1)*[-1,1],[-1,0],m(1)+1,m(1));
D1(1,1) = 0; D1(end,end)=0; D1(end,end-1)=0;

D2 = spdiags(e(2)*[-1,1],[-1,0],m(2)+1,m(2));
D2(1,1) = 0; D2(end,end)=0; D2(end,end-1)=0;

%get normalized flow q and normalized source-field fmat
si   =@(i) prod(h)/h(i);
q{1} = qmat{1}./si(1);
q{2} = qmat{2}./si(2);
q{3} = qmat{3}./si(3);
fmat = Fmat./hd;

%setup cmat
cmat = bsxfun(@times,Cmat,phimat);

%setup divqc (rember that div(q)=fmat
divqc = bsxfun(@times,cmat,fmat);

%setup qGRADc
cTmp  = reshape(cmat,[],k);

m1  = size(q{1});
m2  = size(q{2});
D1c = reshape(D1*cTmp,[m1,k]);
D2c = reshape(D2*cTmp,[m2,k]);

%multiply GRADc and q
qGRADc1 = bsxfun(@times,D1c,q{1});
qGRADc2 = bsxfun(@times,D2c,q{2});

