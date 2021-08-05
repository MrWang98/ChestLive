function [bounds,th]=OtsuMethod(L,ind)

% =========================================================================
% function [bounds,th]=OstuMethod(L,ind)
%
% This function classifies the set of minima curve lengths stored in L into
% two classes by using Otsu's method. Then it returns the boundaries which 
% are supposed to be the meaningful ones.
%
% Inputs:
%   L: set of minima curve lengths
%   ind: original index of each minima
%
% Outputs:
%   bounds: detected bounds
%   th: detected scale threshold
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
% =========================================================================

histo=hist(L,max(L));
%figure;subplot(2,1,1);plot(L);subplot(2,1,2);plot(histo)
Nt=sum(histo);
histo=histo/Nt;

muT=0.0;
for i=1:length(histo)
    muT=muT+i*histo(i);
end

sigbcv=zeros(length(histo)-1,1);

for k=1:(length(histo)-1)
   wb=0.0;
   mu=0.0;
   for i=1:k
       wb=wb+histo(i);
       mu=mu+i*histo(i);
   end
   wf=1-wb;
   mub=mu/wb;
   muf=(muT-mu)/(1-wb);
   sigbcv(k)=wb*wf*(mub-muf)^2;
end

th=maxcheckplateau(sigbcv);

%tag and extract the consistent minima
for i=1:length(L)
   if L(i)<th
       L(i)=0;
   else
       L(i)=1;
   end
end
bounds=ind(find(L));