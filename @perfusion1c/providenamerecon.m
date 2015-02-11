function [name] = providenamerecon(phiopt,Kopt,dim,aiftype,T,deconvmethod,red)

str = perfusion1c.providenameindicator(phiopt,Kopt,dim,aiftype,T);
name = [str '-dcmethod-' deconvmethod '-red-' num2str(red)];

