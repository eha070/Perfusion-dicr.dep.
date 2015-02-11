% CONSTBORDER Placing a constant value around the border of image
%
%   CONSTBORDER(IM,P,V) Placing value V around xy-border as far as P in 
%   image IM. Image must be 2D image, if more than three planes use 
%   constborder3D.
%
function [im] = constborder(im,p,v)

[M N O] = size(im);
blockones = zeros(M,N,O);
blockones(p+1:end-p,p+1:end-p,:) = 1;
blockones = imcomplement(blockones);
im(blockones == 1) = v;
