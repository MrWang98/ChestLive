function [theta,radius] = CreateAngleGrid(W,H)

%=======================================================================
% function [theta,radius] = CreateAngleGrid(W,H)
% 
% This function creates matrices containing the following information for
% each pixel (assuming that the position (0,0) is at the center of the
% image.
%
% Inputs:
%   -W: image width
%   -H: image height
%
% Outputs:
%   -theta : angles in [-pi,pi]
%   -radius: distance from the center
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics and Statistics
% Version: 1.0 - 2013
% Version: 2.0 - 2015
%=======================================================================

theta=zeros(H,W);
radius=zeros(H,W);

middleh=round((H+1)/2);
middlew=round((W+1)/2);

for i=1:W
    for j=1:H
        ri=(i-middlew)*pi/middlew;
        rj=(j-middleh)*pi/middleh;
        radius(j,i)=sqrt(ri^2+rj^2);
        if (ri==0) && (rj==0)
            theta(j,i)=0;
        else
            rt=rj/ri;
            theta(j,i)=atan(rt);
        end
        if ri<0
           if rj<=0
               theta(j,i)=theta(j,i)-pi;
           else
               theta(j,i)=theta(j,i)+pi;
           end
           if theta(j,i)<-3*pi/4
               theta(j,i)=theta(j,i)+2*pi;
           end
        end
    end
end
