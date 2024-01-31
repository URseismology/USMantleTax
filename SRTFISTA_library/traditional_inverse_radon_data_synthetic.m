function [mk1] = traditional_inverse_radon_data_synthetic(d, dt, rayP, q, ...
    flow, fhigh, mu)
%
%  IN   d:         seismic traces
%       dt:        sampling in sec
%       rayP:      ray parameters
%       dist:      epicentral distances
%       q:         q (curvature) range as a vector
%       flow:      freq.  where the inversion starts in Hz (> 0 Hz)
%       fhigh:     freq.  where the inversion ends in Hz (< Nyquist)
%       mu:        regularization parameter
%       alpha:     define alpha
%       tstep:     define tstep
%       maxiter:   maximum iteration count
%       kstop:     steps before maxiter to stop to prevent overfitting
%       vmodel:    [H; Vp; Vs] for travel time curves and t-q points on radon
%       migration: use migration or not when stacking RF
%       twin:      time window in radon and RF plots
%       qwin:      q (curvature) window in radon plots
%
%  OUT  mk1:    sparse radon matrix


% Initialize some parameters
N = 2; % always use parabolic

[nt,np] = size(d);
nq = max(size(q));

nfft = 2*(2^nextpow2(nt));

% initialize model matrices in frequncy domain
D = fft(d,nfft,1);
M0 = zeros(nfft,nq);
i = sqrt(-1);

% frequency set up for kernel initialization
ilow  = floor(flow * dt * nfft) + 1;
if ilow < 2; ilow = 2; end
ihigh = floor(fhigh * dt * nfft) + 1;
if ihigh > floor(nfft/2)+1; ihigh=floor(nfft/2)+1; end
ilow = max(ilow,2);

% initialize kernel matrix (3D) in freq. domain
LLL = zeros(nfft, np, nq);

for ifreq=1:nfft
    f = 2.*pi*(ifreq-1)/nfft/dt;
    LLL(ifreq,:,:) = exp(i*f*(rayP.^N)'*q);
end

% least square initialize model in time domain
Q = eye(nq);

for ifreq=ilow:ihigh
    
    L = squeeze( LLL(ifreq, :,:) );
    y = D(ifreq,:)';
    
    xa = L' * y;
    Ab = L' * L;
    
    A = Ab + mu * Q;
        
    x = A \ xa;
    
    M0(ifreq,:) = x';
    M0(nfft+2-ifreq,:) = conj(x)';
    
end

M0(nfft/2+1,:) = zeros(1,nq);
m0 = real(ifft(M0,[],1));
m0 = m0(1:nt,:);


% % initialize errors
% Dmkm1 = abs(Mkm1); 
mk1 = m0;

return

