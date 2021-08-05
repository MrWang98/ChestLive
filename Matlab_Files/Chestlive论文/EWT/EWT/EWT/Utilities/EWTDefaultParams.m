function params = EWTDefaultParams()

% Perform the detection on the log spectrum instead the spectrum
params.log=1;

params.SamplingRate = -1; %put -1 if you don't know the sampling rate

% Choose the wanted preprocessing for the scales (none,plaw,poly,morpho,tophat)
params.globtrend = 'none';
params.degree=5; % degree for the polynomial interpolation

% Choose the wanted regularization (none,gaussian,average,closing)
params.reg = 'none';
params.lengthFilter = 10;
params.sigmaFilter = 1.5;

% Choose the wanted detection method for the scales (locmax,locmaxmin,ftc,scalespace)
params.detect = 'scalespace';
params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans

params.N = 4; % maximum number of band for the locmaxmin method
params.completion=0; % complete (1) or not (0) the number of boundaries to reach the wanted number.


%% CURVELET PARAMETERS
% Type of curvelet transform (1=scale and radius independent, 2=angles per
% scales)
params.option=1;

% Choose the wanted regularization for the angles
% (none,gaussian,average,closing)
params.curvreg = 'none';
params.curvlengthFilter = 10;
params.curvsigmaFilter = 1.5;

% Choose the wanted global trend removing for the angles (none,plaw,poly,morpho,tophat)
params.curvpreproc='none';
params.curvdegree=4;

% Choose the wanted detection method for the angles (locmax,locmaxmin)
params.curvmethod='scalespace';
params.curvN=2;

