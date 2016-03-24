function [perfmat,lenmat] = flux2perfCMS(qmat,m,h,varargin)
%function CBF = flux2perfCMS(qmat,m,h,varargin)
% CAREFUL, WORKS ONLY FOR m(1) x m(2) x 1-arrays
% working in the right coordinate system was quite some work
%
%INPUT:
% qmat - staggered flow-field (absolute values)
%    m - images-size
%    h - voxel-size
%
%OUTPUT:
% perfmat - perfusion in ml/s/ml
%  lenmat - streamline length in mm
%
%
%
%
%                                      (c)Constantin Sandmann, 18-Mar-2016 
%                                                http://mic.uni-luebeck.de
%                                     See LICENSE.txt for copyright issues
%------------------------------------------------------------------------- 


%default parameters
dt         = min(h)./2;
showFigure = true;

%overwrites default parameter 
for j=1:2:length(varargin),  
	eval([varargin{j},'=varargin{',int2str(j+1),'};']); 
end 

%setup nodal and normalized qx and qy
m       = m(1:2);
[A1,A2] = getInterMatrix(m);
qx      = reshape(A2*qmat{2}(:),m+1)./(h(1)*h(3));
qy      = reshape(-A1*qmat{1}(:),m+1)./(h(2)*h(3));

%setup domain for nodal grid
ygv = fliplr((0:h(1):m(1)*h(1)));
xgv = (0:h(2):m(2)*h(2));
[XNodal,YNodal] = meshgrid(xgv,ygv);

%get domain for cc grid
ygv = fliplr((h(1)/2:h(1):m(1)*h(1)));
xgv = (h(2)/2:h(2):m(2)*h(2));
[Xcc,Ycc] = meshgrid(xgv,ygv);

%set domain
xmin = min(0);
ymin = min(0);
xmax = max(m(2)*h(2));
ymax = max(m(1)*h(1));

%initialize variables
lenmat  = zeros(m); %streamline length
vmat    = zeros(m); %volume of tubes
perfmat = zeros(m); %perfusion

%setup the flow
for i = 1:m(2)
    for j = 1:m(1)
        for k = [1,-1];
            
            %starting position
            posOld  = [Xcc(j,i),Ycc(j,i)];
            POS     = posOld;
            STOP    = false;
            
            %initialize vector with length
            len = [];
            r   = [];                        
            
            %tracking
            while true

                %get direction vector
                vect(1) = interp2(XNodal,YNodal,qx,posOld(1),posOld(2));
                vect(2) = interp2(XNodal,YNodal,qy,posOld(1),posOld(2));
                
                %normalize
                vectn = vect/norm(vect);

                %update position
                pos      = posOld + k*dt*vectn;
                
               %check if new position is on boundary
               invalid = pos(1)<xmin || pos(1)>xmax || ...
                         pos(2)<ymin || pos(2)>ymax;
                      
                %if so: get final update
                if invalid
                    tmp = k*[(xmin-posOld(1))/vectn(1);...
                             (xmax-posOld(1))/vectn(1);...
                             (ymin-posOld(2))/vectn(2);...
                             (ymax-posOld(2))/vectn(2)];
                    dtF = min(tmp(tmp>=0));
                    pos = posOld + k*dtF*vect;
                end
                
               %stopping criterium if arrived at the sink
               if k==1
                   STOP = pos(1) >= max(Xcc(:)) && pos(2) < min(Ycc(:));
               elseif k==-1
                   STOP = pos(1) < min(Xcc(:)) && pos(2) >= max(Ycc(:));
               end               
                
                
                %update length and volume
                r   = [r,1/sqrt(norm(vect))];
                len = [len,norm(pos-posOld)];
                POS = [POS;pos];
                

                if STOP
                    break;
                end                
                
                %update position
                posOld = pos;
                

                
            end;
            
            %write length and volume
            lenmat(j,i)  = lenmat(j,i) + sum(len);
            vmat(j,i)    = vmat(j,i)   + pi*sum((r.^2.*len));   
            perfmat(j,i) = pi./vmat(j,i);
        end 
        
        %show streamline and progress
        if showFigure && any(j==[16;32;48;64])
            figure(1);
            set(1,'name','Perfusion [ml/min/100ml] and parts of streamline');
            colormap parula;
            imagesc(xgv,ygv,perfmat*100*60);
            axis image;
            hold on;
            plot(POS(:,1),POS(:,2),'r','linewidth',2);
            hold off;
            drawnow;
        end

    end;                
end
            
end





%% helper functions

function [A1,A2] = getInterMatrix(m)

    %prepare main matrices
    id  =@(i) speye(m(i)+1);
    %get inter matrices equal bc
    I1          = spdiags(ones(m(1)+1,2),[-1,0],m(1)+1,m(1))./2;
    I1(1,1)     = 1;
    I1(end,end) = 1;
    
    I2          = spdiags(ones(m(2)+1,2),[-1,0],m(2)+1,m(2))./2;
    I2(1,1)     = 1;
    I2(end,end) = 1;
    
    
    %setup avarging matrices
    A2 = kron(id(2),I1);
    A1 = kron(I2,id(1));


end