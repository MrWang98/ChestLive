function Hilb = EWT_InstantaneousComponents(ewt,boundaries)
%==========================================================================
% function Hilb = EWT_InstantaneousComponents(ewt,boundaries)
%
% This function extracts the instantaneous amplitude and phase of each EWT
% component by using the Hilbert transform.
% In order to garantee to have positive frequencies without outliers (e.g.
% the instantaneous phase must be a monotone increasing function), we clean
% the obtained frequencies.
%
% Inputs:
%   -ewt: component obtained by the EWT transform
%   -boundaries: the corresponding boundaries provided by the EWT transform
%
% Outputs:
%   -Hilb: structure containing the instantaneous components (Hilb{i}{1} is
%   the instantaneous amplitude of the component i and Hilb{i}{2} is the
%   instantaneous frequency of component i)
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2014
% Version: 1.0
%==========================================================================



%% compute gamma to estimate the support sizes
gamma=pi;
boundaries=[0 ; boundaries ; pi];
for k=1:length(boundaries)-1
    r=(boundaries(k+1)-boundaries(k))/(boundaries(k+1)+boundaries(k));
    if r<gamma 
       gamma=r;
    end
end

%% We compute the instantaneous amplitude and frequencies for each component
Hilb=cell(length(ewt),2);
for i=1:length(ewt)
    ht=hilbert(ewt{i});

    % Instantaneous amplitude
    Hilb{i}{1}=abs(ht);
   
    %Phase unwrapping, take its derivative and clean the outliers
    Hilb{i}{2}=IFcleaning(diff(unwrap(angle(ht))),(1-gamma)*boundaries(i),(1+gamma)*boundaries(i+1));
       
%     figure;
%     subplot(211);plot(Hilb{i}{1});
%     subplot(212);plot(Hilb{i}{2});
end