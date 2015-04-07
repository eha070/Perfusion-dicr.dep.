clear
clc;
close all;

%load settings
prm = settings;

%create indicatorconv
perfusion1c.createindicatorconv(prm);

%reduce timesteps
aiftype = {'gamma'};
stepred = [25,100,500];
basenameindicator = perfusion1c.providenameindicator(prm.phiopt,prm.Kopt,prm.dim,prm.aiftype,prm.T);
pathload = ['results/syntconv-createindicatorconv' '-' basenameindicator '.mat'];                

for j = 1:numel(stepred);
    perfusion1c.reducetimesteps(pathload,stepred(j));        
end