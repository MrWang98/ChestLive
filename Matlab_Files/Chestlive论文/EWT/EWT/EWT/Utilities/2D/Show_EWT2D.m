function Show_EWT2D(ewt)

%===========================================================
%
% function Show_EWT2D(ewt)
%
% This function permits to plot the outputs of the 2D EWT.
%
% Input:
%   - ewt: cell containing the EWT outputs
%
% Author: J.Gilles
% Institution: UCLA - Department of Mathematics
% email: jegilles@math.ucla.edu
% Date: March, 1st, 2013
%
%===========================================================
p=ceil(length(ewt)/2);

figure;
for n=1:length(ewt)
    subplot(2,p,n);imshow(ewt{n},[]);
end