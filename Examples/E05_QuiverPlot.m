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

load('synt-createflowTPFA-phi-flat-K-flat-dim-64.mat')


%%
%setup domain
omega = [0,1,0,1];
m     = size(pmat);


%% convert flow

%get flow
qcc = perfusion1c.convertFlowStagToCC(qmat);
qy  = qcc{1};
qx  = qcc{2};



%% show 2D-streamlines

%get meshgrid to plot qcc, the cell-centered flow
ygv   = linspace(omega(2),omega(1),m(1)); %from (0 to 1)
xgv   = linspace(omega(3),omega(4),m(2)); %this goes from (1 to 0) (matlab, duh)
[x,y] = meshgrid(xgv,ygv);



%setup points where to follow the streamlines
startx = (0:.1:1);
starty = .9*ones(size(startx));



%create quiver plot. Notice that flow in y-direction is reversed since the
%grid is reversed
figure(1);clf;
quiver(x(:),y(:),qx(:),-qy(:),10)
streamline(x,y,qx,-qy,startx,starty)
xlabel('x');
ylabel('y');
set(gca,'FontSize',15);


%% show 3D-streamlines

%number of timepoints to simulate
k = 100;

%start-coordinates for the streamlines
xstart = .1;
ystart = .9;
tstart = 0;
numt   = 1;

tend   = 10;



%setup starting points for streamlines
startx = xstart*ones(numt,1);
starty = ystart*ones(numt,1);
startz = linspace(tstart,tend,numt)';

%setup stepsizes to downsample the data for visualization
stepxy = 5;
stepz = 20;

%get 3D-flow
qx3D = repmat(qx,[1,1,k]);
qy3D = repmat(qy,[1,1,k]);
qz3D = ones([m,k]);

%get meshgrid to plot qcc, the cell-centered flow
ygv   = linspace(omega(2),omega(1),m(1)); %from (0 to 1)
xgv   = linspace(omega(3),omega(4),m(2)); %this goes from (1 to 0) (matlab, duh)
zgv   = linspace(0,60,k);
[x,y,z] = meshgrid(xgv,ygv,zgv);


%downsample the vector-fields for display
idxX = (1:stepxy:m(2));
idxY = (1:stepxy:m(1));
idxZ = (1:stepz:k);

qxLow = qx3D(idxY,idxX,idxZ);
qyLow = qx3D(idxY,idxX,idxZ);
qzLow = qx3D(idxY,idxX,idxZ);
xLow  = x(idxY,idxX,idxZ);
yLow  = y(idxY,idxX,idxZ);
zLow  = z(idxY,idxX,idxZ);



%create quiver plot with streamlines. Notice that flow in y-direction is reversed since the
%grid is reversed
figure(2);clf;
hq = quiver3(xLow,yLow,zLow,qxLow,-qyLow,qzLow,1.5);
hs = streamline(x,y,z,qx3D,-qy3D,qz3D,startx,starty,startz);
xlabel('x');
ylabel('y');
zlabel('time');
set(hs,'LineWidth',2);
set(gca,'FontSize',15);


