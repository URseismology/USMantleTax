function [mk1, errors] = sparse_inverse_radon_data_synthetic(d, dt, rayP, q, ...
    flow, fhigh, mu, alpha, tstep, maxiter, kstop, gt)
%SPARSE_INVERSE_RADON_DATA_SYNTHETIC Calculates SRTIS-based reconstruction
%on synthetic data.
%   Inputs:
%       d:         seismic traces
%       dt:        sampling in sec
%       rayP:      ray parameters
%       q:         q (curvature) range as a vector
%       flow:      freq.  where the inversion starts in Hz (> 0 Hz)
%       fhigh:     freq.  where the inversion ends in Hz (< Nyquist)
%       mu:        regularization parameter of the starting point
%       alpha:     define alpha
%       tstep:     define tstep
%       maxiter:   maximum iteration count
%       kstop:     steps before maxiter to stop to prevent overfitting
%       gt:        true Radon image
%
%   Outputs:
%       mk1:       sparse Radon image
%       errors:    error between the output of the algorithm and Radon
%           image at every iteration, i.e., convergence curve
%

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
LLLINV = zeros(nfft, nq, nq);

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
    
    LLLINV(ifreq,:,:) =  pinv(Ab + 1.0 * Q);
    
    x = A \ xa;
    
    M0(ifreq,:) = x';
    M0(nfft+2-ifreq,:) = conj(x)';
    
end

M0(nfft/2+1,:) = zeros(1,nq);
m0 = real(ifft(M0,[],1));
m0 = m0(1:nt,:);

% setup initial models
Mkm1 = M0; % mk-1 , k = 1
mkm1 = m0; % mk-1 , k = 1

% % initialize errors
% Dmkm1 = abs(Mkm1); 
mk1 = m0;
% dkm1 = zeros(size(d));
% 
% taus = (0:1:nt-1) * dt; %generate taus scaled by the sample time

% % parse velocity model
% H  = vmodel(1, :);
% Vp = vmodel(2, :);
% Vs = vmodel(3, :);

% % calculate arrivals
% shft = 1;
% [t1, t2, t3] = travelTimesAppx(Vp, Vs, H, rayP, 1, 2);
% Tphase = [t1, t2, t3] - shft;

% % calculate tau-q points
% [all_q, all_tau] = get_q_t(H, Vp, Vs);
% tqpts_posq = [all_tau(:,1) all_q(:,1)]; % all direct
% tqpts_negq = [[all_tau(:,2); all_tau(:,3)] [all_q(:,2); all_q(:,3)]]; % all multiples

% % get diagonal mask
% [Kmask, qps_q] = KdiagMask(m0', taus, q, H, Vp, Vs, 2);
% Kmask = Kmask';

errors = zeros(1, maxiter+1);
errors(1) = sum(abs(mk1-gt),"all");

% start iteration
for kstep = 1:maxiter
    
%     % if last step, mask radon image and do forward
%     if kstep == maxiter - kstop
%         [Mk1, ~] = shrink_tdsrt_mask(mkm1, Dmkm1, kstep, ...
%             maxiter, alpha,tstep, nfft, Kmask); % add in the update after shrinking and fourier ...
%         [~, dkm1] = radon3d_forward(LLL, Mk1, nt, ilow, ihigh); % transform LSRT model into data (D) domain -- freq.
%     end
%     
%     rawRF = d';
%     filteredRF = dkm1';
%     
%     if migration
%         [~, rawRF] = migrateAndStackRF(rawRF, tt, rayP, H, Vp, Vs, target);
%         [~, filteredRF] = migrateAndStackRF(filteredRF, tt, rayP, H, Vp, Vs, target);
%     end
%     
%     rawRF_stack = mean(rawRF, 1);
%     rawRF_stack = rawRF_stack / max(abs(rawRF_stack));
%     filteredRF_stack = pws_stack(rawRF, filteredRF);
%     
%     
%     %%%%%%%%%%%%%%%%%%
%     % start plotting %
%     %%%%%%%%%%%%%%%%%%
%     figure(15);
%     clf;
%     set(gcf,'position',[50,50,1200,1200]);
%     
%     %%%%% top left: sparse radon %%%%%
%     subplot(10,10,[1:4,11:14,21:24,31:34]);
%     % plot radon image
%     RadonPlot(taus, q, mk1, twin, qwin, tqpts_posq - shft, tqpts_negq - shft);
%     % plot mask before apply
%     hold on;
%     plot(taus, qps_q, 'r')
%     plot(taus, qps_q - 50, 'r--')
%     plot(taus, qps_q + 50, 'r--')
%     
%     ylabel('$q (s/km)^2$', 'interpreter', 'latex', 'Fontsize', 20 );
%     title(['$m_{ ' num2str(kstep) '} = \mathcal{R^{\dag}}_{sp}(d)$'], 'interpreter', 'latex', 'fontsize', 20);
%     
%     %%%%% middle left: filtered RF %%%%%
%     subplot(10,10,[41:44,51:54,61:64,71:74,81:84]);
%     RFWigglePlot(filteredRF, taus, rayP, dist, Tphase, twin, 0.5, 0,1);
%     if kstep == maxiter - kstop
%         title(['$\hat{d}_{ ', num2str(kstep) '} = \mathcal{R}(m_{' num2str(kstep) '}K_{q}^{+})$'], 'interpreter', 'latex', 'fontsize', 20);
%     else
%         title(['$\hat{d}_{ ', num2str(kstep) '} = \mathcal{R}(m_{' num2str(kstep) '})$'], 'interpreter', 'latex', 'fontsize', 20);
%     end
%     
%     %%%%% middle right: raw RF %%%%%
%     subplot(10,10,[47:50,57:60,67:70,77:80,87:90])
%     RFWigglePlot(rawRF, taus, rayP, dist, Tphase, twin, 0.5, 0,1);
%     title('$d$', 'interpreter', 'latex', 'fontsize', 20);
%     
%     %%%%% top right: initial radon from LS %%%%%
%     subplot(10,10,[7:10,17:20,27:30,37:40]);
%     RadonPlot(taus, q, m0, twin, qwin, [], []);
%     title('$m_{0} =\mathcal{R^{\dag}}_{ls}(d)$',  'interpreter', 'latex', 'Fontsize', 20);
%     
%     %%%%% bottom left: stack filtered RF %%%%%
%     subplot(10,10,91:94);
%     jbfill(taus, max(filteredRF_stack, 0), zeros(1, length(taus)), [0 0 1],'k', 1, 1.0);
%     jbfill(taus, min(filteredRF_stack, 0), zeros(1, length(taus)), [1 0 0],'k', 1, 1.0);
%     xlim(twin);
%     ylim([-1 1]);
%     xlabel('Time (s)','FontSize', 20);
%     
%     %%%%% bottom right: stack raw RF %%%%%
%     subplot(10,10,97:100);
%     jbfill(taus, max(rawRF_stack, 0), zeros(1, length(taus)), [0 0 1],'k', 1, 1.0);
%     jbfill(taus, min(rawRF_stack, 0), zeros(1, length(taus)), [1 0 0],'k', 1, 1.0);
%     xlim(twin);
%     ylim([-1 1]);
%     xlabel('Time (s)','FontSize', 20)
    
%     if kstep == maxiter - kstop
%         break;
%     end
    
    % do forward and invserse
    [Dkm1, dkm1] = radon3d_forward(LLL, Mkm1, nt, ilow, ihigh); % transform LSRT model into data (D) domain -- freq.
    DDkm1 = D - Dkm1;       % prediction error in data domain -- Freq
    
    [~, Dmkm1] = radon3d_pinv(LLL, LLLINV, DDkm1, nt, ilow, ihigh); % function to do pseudoinv(L'L)L'*[x] in 3D
    
    [Mk1, mk1] = shrink_tdsrt(mkm1, Dmkm1, kstep, ...
        maxiter, alpha,tstep, nfft); % add in the update after shrinking and fourier ...
    
    Mkm1 = Mk1;  % step up
    mkm1 = mk1;  % step up and loop

    errors(kstep+1) = sum(abs(mk1-gt),"all");
    
    pause(.0002)
    
end

return

