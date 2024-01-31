function [time, binvals, radRF, transRF, varR, varT] = ...
    MTCRF_ET(RR, TT, ZZ, NN, Dt, Fcut, binparams, dist, rayp, taperWin, olap)


% get binning parameters
binmin = binparams(1);
binmax = binparams(2);
bindelta = binparams(3);

hdist = histogram(dist, 'BinMethod', 'sqrt');

stackparams.binedge = binmin:bindelta:binmax;
stackparams.nbins = length(stackparams.binedge)-1;
stackparams.deltabin = binparams(4);

stackparams.bincnt = hdist.Values;
stackparams.binlims = hdist.BinLimits;
stackparams.distvals = hdist.Data;

% define ray parameter (x) and epicentral distance (y) column vectors for
% the stacked receiver functions
x = dist;
if isrow(x)
    x = x';
end

y = rayp;
if isrow(y)
    y = y';
end

% find unique values of ray parameter (x) and corresponding epicentral
% distance (y)
[xu, iu, ~] = unique(x);
yu = y(iu);

% define x and y bins as column vectors
xbin = stackparams.binedge(1:stackparams.nbins);
ybin = interp1(xu, yu, xbin, 'linear', 'extrap');
binvals = [xbin; ybin]';

% initialize output matrices for RF and variance
[N, nwaves] = size(RR);
radRF = zeros(N, nwaves);
transRF = zeros(N, nwaves);
varR = zeros(N, nwaves);
varT = zeros(N, nwaves);


% pre-define slepian taper using input taper window length
P = 2.5;
K = 2;
Ntaper = round(taperWin / Dt);
[PSI, ~] = sleptap(Ntaper, P, K);

% for loop to calculate single trace RF

for iwave = 1: nwaves
    
    R = RR(:, iwave);
    T = TT(:, iwave);
    Z = ZZ(:, iwave);
    dnoi = NN(:, iwave);
    
    [~, radRF(:, iwave), transRF(:,iwave), varR(:, iwave), varT(:,iwave)] = ...
        MTCRF_ET_se(R, T, Z, dnoi, PSI, K, Dt, taperWin, olap);
    
end

% stack
[radRFstck, transRFstck] = stackwithovlpRFs(radRF, transRF, varR, varT, stackparams);
[ntime, nbins] = size(radRFstck);

% time domain receiver functions after stackign
radRF = zeros(ntime, nbins);
transRF = zeros(ntime, nbins);

% transform to time domain
fmax = 1/(2.0*Dt);

for ibin = 1: nbins
    
    HRin = radRFstck(:, ibin);
    HTin = transRFstck(:, ibin);
    
    [HR, ~] = costaper(HRin, Fcut, Dt);
    [HT, ~] = costaper(HTin, Fcut, Dt);
    
    radRF(:, ibin) = (2*fmax/Fcut) .* fftshift(real(ifft(HR)));  % -tlag to +tlag in time
    transRF(:, ibin) = (2*fmax/Fcut) .* fftshift(real(ifft(HT)));
    
end

time = ((0:N-1)-floor(N/2)) * Dt;
