function noisyRF = RFaddNoise(RFsyn,tt,s2n,localBaseDir)

% adding realistic noise to a recieer function
addpath([localBaseDir '/softwares/DOTM/']);  % bandpass
addpath([localBaseDir '/softwares/11_crewes/syntraces/']); % rnoise

nsample = length(tt);

dt = diff(tt); sampt = dt(1); Fs = 1./sampt; %signal sampling 

nse = rnoise(RFsyn,s2n, 1:nsample, 0);
nsefilt = bandpass(nse, Fs, 0.01, 0.2, 2);

pn1= sqrt(nse*nse');
pn2= sqrt(nsefilt'*nsefilt);
scalar = pn1/(pn2*1);            % data variance preserced

nsefilt = nsefilt .* scalar;
noisyRF = RFsyn + nsefilt';    % add to RF here

rmpath([localBaseDir '/softwares/DOTM/']);  % bandpass
rmpath([localBaseDir '/softwares/11_crewes/syntraces/']); % rnoise
