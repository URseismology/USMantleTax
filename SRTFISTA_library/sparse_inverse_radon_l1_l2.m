function m = sparse_inverse_radon_l1_l2(d, dt, rayP, q, ...
    flow, fhigh, mu, alpha, lambd, rho, maxiter)
%%SPARSE_INVERSE_RADON_L1_L2 Calculates SRTL1-2-based reconstruction on
%%real data
%
% Inputs:
%       d:         seismic traces
%       dt:        sampling in sec
%       rayP:      ray parameters
%       q:         q (curvature) range as a vector
%       flow:      freq.  where the inversion starts in Hz (> 0 Hz)
%       fhigh:     freq.  where the inversion ends in Hz (< Nyquist)
%       mu:        regularization parameter of the L2 regularized solution
%            used as the starting point of the algorithm
%       alpha:     regularization parameter in front of the L1
%            regularization
%       lambda:    regularization parameter in front of the L1-L2 term
%       rho:       penalty parameter of the augmented Lagrangian
%       maxiter:   maximum number of iterations 
%
% Outputs:
%       m:         reconstructed Radon image 
%


%% Initialization
% I copied this part from Evan's code

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


%% Perform projection on the data
% In Evan's code, L2-regularized solution is calculated, which requires
% tuning another parameter. To remove the parameter, I have decided to
% perform "poor man's inversion", which is A^H*d. It seems to work fine as
% a starting point.

LLLINV = zeros(nfft, nq, nq);
% least square initialize model in time domain
Q = eye(nq);
for ifreq=ilow:ihigh
    L = squeeze( LLL(ifreq, :,:) );
    y = D(ifreq,:)';

    xa = L' * y;
    Ab = L' * L;

    A = Ab + mu * Q;
    LLLINV(ifreq,:,:) = (Ab + rho * Q);

    x = A \ xa;

    M0(ifreq,:) = x';
    M0(nfft+2-ifreq,:) = conj(x)';

end
M0(nfft/2+1,:) = zeros(1,nq);
m0 = real(ifft(M0,[],1));
m0 = m0(1:nt,:);


%% Some intermediate calculations
% I copied this part from Evan's code.

% taus = (0:1:nt-1) * dt; %generate taus scaled by the sample time


%% L1-L2 algorithm
% Main algorithm. Note that I removed a lot of the lines used for plotting
% purposes since my data is simulated and I do not have some of the
% required variables.

% initialize the algorithm.
m = m0;
z = m;
w = ones(size(z));


for kstep = 1:maxiter

%     if verbose == 1
%         %%%%%%%%%%%%%%%%%%
%         % start plotting %
%         %%%%%%%%%%%%%%%%%%
%         figure(15);
%         clf;
%         set(gcf,'position',[50,50,1200,1200]);
% 
%         %%%%% top left: sparse radon %%%%%
%         subplot(10,10,[1:4,11:14,21:24,31:34]);
%         % plot radon image
%         RadonPlot(taus, q, m, twin, qwin, [], []);
%         % plot mask before apply
%         hold on;
% 
%         ylabel('$q (s/km)^2$', 'interpreter', 'latex', 'Fontsize', 20 );
%         title(['$m_{ ' num2str(kstep) '} = \mathcal{R^{\dag}}_{sp}(d)$'], 'interpreter', 'latex', 'fontsize', 20);
% 
%         %%%%% top right: initial radon from LS %%%%%
%         subplot(10,10,[7:10,17:20,27:30,37:40]);
%         RadonPlot(taus, q, m0, twin, qwin, [], []);
%         title('$m_{0} =\mathcal{R^{\dag}}_{ls}(d)$',  'interpreter', 'latex', 'Fontsize', 20);
% 
%     end

    % calculate F[g_k]
    Fg = alpha * lambd * ifft(real(m),nfft,1) / sqrt(sum(abs(fft(real(m),nfft,1)).^2, "all"));
    % m-update
    temp1 = - Fg + rho * fft(real(z),nfft,1) - fft(real(w),nfft,1);
    [~, m] = radon3d_pinv_2(LLL, LLLINV, fft(real(d),nfft,1), temp1, nt, ilow, ihigh);
    % z-update
    z = wthresh(m + w ./ rho,'s',lambd / rho);
    % w-update
    w = w + m - z;


    pause(.0002)

end










