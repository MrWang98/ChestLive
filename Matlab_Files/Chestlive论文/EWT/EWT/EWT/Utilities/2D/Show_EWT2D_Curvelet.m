function Show_EWT2D_Curvelet(ewtc)

%==========================================================================
% function Show_EWT2D_Curvelet(ewtc)
% 
% This function displays the curvelet coefficient obtained by the Empirical 
% Curvelet Transform (one figure per scale)
%
% Inputs:
%   -ewtc: output of the empirical curvelet transform
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2013
% Version: 1.0
%==========================================================================

Ns=length(ewtc);

% Show the low pass filter
figure;
imshow(ewtc{1},[]);

% mincur=min(ewtc{2}{1}(:));
% maxcur=max(ewtc{2}{1}(:));
% 
% for s=2:Ns
%    Nt=length(ewtc{s});
%    for t=1:Nt
%         mincur=min(mincur,min(ewtc{s}{t}(:)));
%         maxcur=max(maxcur,max(ewtc{s}{t}(:)));
%    end
% end
% maxcur=maxcur-mincur;

for s=2:Ns
   Nt=length(ewtc{s});
   figure;
   Nr=ceil(sqrt(Nt));
   for t=1:Nt
      subplot(Nr,Nr,t);
%      imshow((ewtc{s}{t}-min(ewtc{s}{t}(:)))/(max(ewtc{s}{t}(:))-min(ewtc{s}{t}(:)))*maxcur-mincur,[]);
      imshow(ewtc{s}{t},[]);
   end
end