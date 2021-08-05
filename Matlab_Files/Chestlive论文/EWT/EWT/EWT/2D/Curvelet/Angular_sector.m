function Angular = Angular_sector(theta,radius,theta0,theta1,r0,r1,gammaw,Dtheta,last)

%==================================================================================
% function Angular = Angular_sector(theta,radius,theta0,theta1,r0,r1,gammaw,Dtheta)
% 
% This function creates the curvelet filter in the Fourier domain which
% have a support defined in polar coordinates (r,angles) in
% [(1-gammaw)r0,(1+gammaw)r1]x[theta0-Dtheta,theta1+Dtheta] and 
% [(1-gammaw)r0,(1+gammaw)r1]x[theta0-Dtheta+pi,theta1+Dtheta+pi]
%
% Inputs:
%   -theta: angle grid
%   -radius: radius grid
%   -theta0,theta1: the consecutive angles defining the support's limits
%   (theta0 < theta1)
%   -r0,r1: the consecutive scales defining the support's limits
%   (r0<r1<=pi)
%   -gammaw: coefficient defining the spread of scale transition zones
%   -Dtheta: coefficient defining the spread of angle transition zones
%   -last: 1 if we consider the last angular sector, 0 otherwise
%
% Output:
%   -Angular: curvelet filter in Fourier domain
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics and Statistics
% Version: 1.0 (2013)
% Version: 2.0 (2015)
%==================================================================================

H=size(theta,2);
W=size(theta,1);

mj=round((W+1)/2);
mi=round((H+1)/2);

wan=1/(2*gammaw*r0);
wam=1/(2*gammaw*r1);
wpbn=(1+gammaw)*r0;
wmbn=(1-gammaw)*r0;
wpbm=(1+gammaw)*r1;
wmbm=(1-gammaw)*r1;

an=1/(2*Dtheta);
pbn=theta0+Dtheta;
mbn=theta0-Dtheta;
pbm=theta1+Dtheta;
mbm=theta1-Dtheta;


Angular=zeros(size(theta));

if r1<pi
   for i=1:H
       for j=1:W
            if ((theta(j,i)>pbn) && (theta(j,i)<mbm))
                if ((radius(j,i)>=wpbn) && (radius(j,i)<=wmbm)) %inside
                    Angular(j,i)=1;
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>wmbm) && (radius(j,i)<=wpbm)) %top of inside - radial only
                    Angular(j,i)=cos(pi*EWT_beta(wam*(radius(j,i)-wmbm))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom of inside - radial only
                    Angular(j,i)=sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                end
            elseif ((theta(j,i)>=mbn) && (theta(j,i)<=pbn)) 
                if ((radius(j,i)>=wpbn) && (radius(j,i)<=wmbm)) %left of inside - angular only
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>wmbm) && (radius(j,i)<=wpbm)) %top-left of inside - radial/angular mix
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2)*cos(pi*EWT_beta(wam*(radius(j,i)-wmbm))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-left of inside - radial/angular mix
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                end
            elseif ((theta(j,i)>=mbm) && (theta(j,i)<=pbm))
                if ((radius(j,i)>=wpbn) && (radius(j,i)<=wmbm)) %right of inside - angular only
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>wmbm) && (radius(j,i)<=wpbm)) %top-right of inside - radial/angular mix
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2)*cos(pi*EWT_beta(wam*(radius(j,i)-wmbm))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-right of inside - radial/angular mix
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                end
            end
       end
   end
else
   if rem(W,2)==0
       ii=2;
   else
       ii=1;
   end
   if rem(H,2)==0
       jj=2;
   else
       jj=1;
   end
    
   for i=ii:H
       for j=jj:W
            if ((theta(j,i)>pbn) && (theta(j,i)<mbm))
                if (radius(j,i)>=wpbn) %inside
                    Angular(j,i)=1;
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom of inside - radial only
                    Angular(j,i)=sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                end
             elseif ((theta(j,i)>=mbn) && (theta(j,i)<=pbn))
                if (radius(j,i)>=wpbn) %left of inside - angular only
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-left of inside - radial/angular mix
                    Angular(j,i)=sin(pi*EWT_beta(an*(theta(j,i)-mbn))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                end
            elseif ((theta(j,i)>=mbm) && (theta(j,i)<=pbm)) 
                if (radius(j,i)>=wpbn) %right of inside - angular only
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                elseif ((radius(j,i)>=wmbn) && (radius(j,i)<=wpbn)) %bottom-right of inside - radial/angular mix
                    Angular(j,i)=cos(pi*EWT_beta(an*(theta(j,i)-mbm))/2)*sin(pi*EWT_beta(wan*(radius(j,i)-wmbn))/2);
                    Angular(2*mj-j,2*mi-i)=Angular(j,i);
                end
            end
       end
   end

    if rem(W,2)==0 %process the first column
        for j=1:H
            the=atan((j-mj)*mi/((1-mi)*mj));
            if ((the>pbn) && (the<mbm))
                if (radius(j,1)>=wpbn) %inside
                    Angular(j,1)=1;
                elseif ((radius(j,1)>=wmbn) && (radius(j,1)<=wpbn)) %bottom of inside - radial only
                    Angular(j,1)=sin(pi*EWT_beta(wan*(radius(j,1)-wmbn))/2);
                end
            elseif ((the>=mbn) && (the<=pbn))
                if (radius(j,1)>=wpbn) %left of inside - angular only
                    Angular(j,1)=sin(pi*EWT_beta(an*(the-mbn))/2);
                elseif ((radius(j,1)>=wmbn) && (radius(j,1)<=wpbn)) %bottom-left of inside - radial/angular mix
                    Angular(j,1)=sin(pi*EWT_beta(an*(the-mbn))/2)*sin(pi*EWT_beta(wan*(radius(j,1)-wmbn))/2);
                end
            elseif ((the>=mbm) && (the<=pbm)) 
                if (radius(j,1)>=wpbn) %right of inside - angular only
                    Angular(j,1)=cos(pi*EWT_beta(an*(the-mbm))/2);
                elseif ((radius(j,1)>=wmbn) && (radius(j,1)<=wpbn)) %bottom-right of inside - radial/angular mix
                    Angular(j,1)=cos(pi*EWT_beta(an*(the-mbm))/2)*sin(pi*EWT_beta(wan*(radius(j,1)-wmbn))/2);
                end
            end   
        end
    end
    
    if rem(H,2)==0 %process the first row
        if last==1
            theta(1,:)=theta(1,:)+pi;
        end
        
        for i=2:W
            if ((theta(1,i)>pbn) && (theta(1,i)<mbm))
                if (radius(1,i)>=wpbn) %inside
                    Angular(1,i)=1;
                elseif ((radius(1,i)>=wmbn) && (radius(1,i)<=wpbn)) %bottom of inside - radial only
                    Angular(1,i)=sin(pi*EWT_beta(wan*(radius(1,i)-wmbn))/2);
                end
            elseif ((theta(1,i)>=mbn) && (theta(1,i)<=pbn))
                if (radius(1,i)>=wpbn) %left of inside - angular only
                    Angular(1,i)=sin(pi*EWT_beta(an*(theta(1,i)-mbn))/2);
                elseif ((radius(1,i)>=wmbn) && (radius(1,i)<=wpbn)) %bottom-left of inside - radial/angular mix
                    Angular(1,i)=sin(pi*EWT_beta(an*(theta(1,i)-mbn))/2)*sin(pi*EWT_beta(wan*(radius(1,i)-wmbn))/2);
                end
            elseif ((theta(1,i)>=mbm) && (theta(1,i)<=pbm)) 
                if (radius(1,i)>=wpbn) %right of inside - angular only
                    Angular(1,i)=cos(pi*EWT_beta(an*(theta(1,i)-mbm))/2);
                elseif ((radius(1,i)>=wmbn) && (radius(1,i)<=wpbn)) %bottom-right of inside - radial/angular mix
                    Angular(1,i)=cos(pi*EWT_beta(an*(theta(1,i)-mbm))/2)*sin(pi*EWT_beta(wan*(radius(1,i)-wmbn))/2);
                end
            end
         end
     end
end