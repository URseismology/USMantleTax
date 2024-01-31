function m = sparse_inverse_radon_fista(d, dt, rayP, q, ...
    flow, fhigh, reg_param, reg_param2, maxiter)
%%SPARSE_INVERSE_RADON_FISTA Calculates SRTFISTA-based reconstruction on
%%real data
%
% Inputs:
%       d:         seismic traces
%       dt:        sampling in sec
%       rayP:      ray parameters
%       q:         q (curvature) range as a vector
%       flow:      freq.  where the inversion starts in Hz (> 0 Hz)
%       fhigh:     freq.  where the inversion ends in Hz (< Nyquist)
%       reg_param: regularization parameter (denoted by lambda in the notes)
%       reg_param2:regularization parameter of the L2 regularized solution
%            used as the starting point of the algorithm
%       maxiter:   maximum number of iterations 
%
% Outputs:
%       m:         reconstructed Radon image 
%


%% Initialization

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
% Calculate the starting point of the algorithm with L2 regularization.

Q = eye(nq);

for ifreq=ilow:ihigh

    L = squeeze( LLL(ifreq, :,:) );
    y = D(ifreq,:)';

    xa = L' * y;
    Ab = L' * L;

    A = Ab + reg_param2 * Q;

    x = A \ xa;

    M0(ifreq,:) = x';
    M0(nfft+2-ifreq,:) = conj(x)';

end

M0(nfft/2+1,:) = zeros(1,nq);
m0 = real(ifft(M0,[],1));
m0 = m0(1:nt,:);


%% Calculate the step size
% To calculate the step size, we estimate the maximum eigenvalue of A^H*A
% using power iterations. Here, L is the resulting estimate. Here, we only
% perform 2 iterations, but it can be increased to improve the accuracy of
% the estimate.
%
% https://en.wikipedia.org/wiki/Power_iteration

b_k = rand(nt, nq);
B_k = fft(real(b_k),nfft,1);
for i=1:2
    [B_k1, ~] = radon3d_forward(LLL, B_k, nt, ilow, ihigh);
    [~, b_k1] = radon3d_forward_adjoint(LLL, B_k1, nt, ilow, ihigh);
    b_k1_norm = sqrt(sum(b_k1(:).^2));
    b_k = b_k1 ./ b_k1_norm;
    B_k = fft(real(b_k),nfft,1);
end
[B_k_temp, ~] = radon3d_forward(LLL, B_k, nt, ilow, ihigh);
[~, b_k_temp] = radon3d_forward_adjoint(LLL, B_k_temp, nt, ilow, ihigh);
L = sum((b_k .* b_k_temp),"all") / sum((b_k .* b_k),"all");


% %% Some intermediate calculations
% 
% taus = (0:1:nt-1) * dt; %generate taus scaled by the sample time
% 
% % parse velocity model
% H  = vmodel(1, :);
% Vp = vmodel(2, :);
% Vs = vmodel(3, :);
% 
% % calculate arrivals
% shft = 1;
% [t1, t2, t3] = travelTimesAppx(Vp, Vs, H, rayP, 1, 2);
% Tphase = [t1, t2, t3] - shft;
% 
% % calculate tau-q points
% [all_q, all_tau] = get_q_t(H, Vp, Vs);
% tqpts_posq = [all_tau(:,1) all_q(:,1)]; % all direct
% tqpts_negq = [[all_tau(:,2); all_tau(:,3)] [all_q(:,2); all_q(:,3)]]; % all multiples
% 
% % get diagonal mask
% [Kmask, qps_q] = KdiagMask(m0', taus, q, H, Vp, Vs, 2);
% Kmask = Kmask';



%% FISTA
% Main algorithm. Note that I removed a lot of the lines used for plotting
% purposes since my data is simulated and I do not have some of the
% required variables.

% initialize the algorithm.
m = m0;
s = m;
q_t = 1;
step_size = 1 / L * 0.9;

% obj_function_vals = [];

for kstep = 1:maxiter

    % % apply the mask in the last iteration
    %if kstep == maxiter
    %    m = m .* Kmask;
    %end

%     % show the plots if needed
%     if verbose == 1
%         [~, dkm1] = radon3d_forward(LLL, fft(real(m),nfft,1), nt, ilow, ihigh); % transform LSRT model into data (D) domain -- freq.
% 
%         rawRF = d';
%         filteredRF = dkm1';
% 
%         if migration
%             [~, rawRF] = migrateAndStackRF(rawRF, tt, rayP, H, Vp, Vs, target);
%             [~, filteredRF] = migrateAndStackRF(filteredRF, tt, rayP, H, Vp, Vs, target);
%         end
% 
%         rawRF_stack = mean(rawRF, 1);
%         rawRF_stack = rawRF_stack / max(abs(rawRF_stack));
%         filteredRF_stack = pws_stack(rawRF, filteredRF);
% 
%         figure(1);
%         clf;
%         set(gcf,'position',[50,50,1200,1200]);
% 
%         %%%%% top left: sparse radon %%%%%
%         subplot(10,10,[1:4,11:14,21:24,31:34]);
%         % plot radon image
%         RadonPlot(taus, q, m, twin, qwin, tqpts_posq - shft, tqpts_negq - shft);
%         % plot mask before apply
%         hold on;
%         plot(taus, qps_q, 'r')
%         plot(taus, qps_q - 50, 'r--')
%         plot(taus, qps_q + 50, 'r--')
% 
%         ylabel('$q (s/km)^2$', 'interpreter', 'latex', 'Fontsize', 20 );
%         title(['$m_{ ' num2str(kstep) '} = \mathcal{R^{\dag}}_{sp}(d)$'], 'interpreter', 'latex', 'fontsize', 20);
% 
%         %%%%% middle left: filtered RF %%%%%
%         subplot(10,10,[41:44,51:54,61:64,71:74,81:84]);
%         RFWigglePlot(filteredRF, taus, rayP, dist, Tphase, twin, 0.5, 0,1);
%         title(['$\hat{d}_{ ', num2str(kstep) '} = \mathcal{R}(m_{' num2str(kstep) '})$'], 'interpreter', 'latex', 'fontsize', 20);
% 
%         %%%%% middle right: raw RF %%%%%
%         subplot(10,10,[47:50,57:60,67:70,77:80,87:90])
%         RFWigglePlot(rawRF, taus, rayP, dist, Tphase, twin, 0.5, 0,1);
%         title('$d$', 'interpreter', 'latex', 'fontsize', 20);
% 
%         %%%%% top right: initial radon from LS %%%%%
%         subplot(10,10,[7:10,17:20,27:30,37:40]);
%         RadonPlot(taus, q, m0, twin, qwin, [], []);
%         title('$m_{0} =\mathcal{R^{\dag}}_{ls}(d)$',  'interpreter', 'latex', 'Fontsize', 20);
% 
%         %%%%% bottom left: stack filtered RF %%%%%
%         subplot(10,10,91:94);
%         jbfill(taus, max(filteredRF_stack, 0), zeros(1, length(taus)), [0 0 1],'k', 1, 1.0);
%         jbfill(taus, min(filteredRF_stack, 0), zeros(1, length(taus)), [1 0 0],'k', 1, 1.0);
%         xlim(twin);
%         ylim([-1 1]);
%         xlabel('Time (s)','FontSize', 20);
% 
%         %%%%% bottom right: stack raw RF %%%%%
%         subplot(10,10,97:100);
%         jbfill(taus, max(rawRF_stack, 0), zeros(1, length(taus)), [0 0 1],'k', 1, 1.0);
%         jbfill(taus, min(rawRF_stack, 0), zeros(1, length(taus)), [1 0 0],'k', 1, 1.0);
%         xlim(twin);
%         ylim([-1 1]);
%         xlabel('Time (s)','FontSize', 20)
% 
%         %%%%% loss curve %%%%
%         [~, temp] = radon3d_forward(LLL, fft(real(m),nfft,1), nt, ilow, ihigh);
%         data_fidelity = sqrt(sum(abs(temp-d).^2, "all"));
%         prior = reg_param * sum(abs(m), "all");
%         obj_func_val = data_fidelity + prior;
%         obj_function_vals = [obj_function_vals, obj_func_val];
% 
%         figure(2);
%         semilogy(1:kstep, obj_function_vals)
%         title('Value of the objective function')
%         ylabel('0.5*||F^H*L*F*m^(k) - d ||_2^2 + lambda*||m^(k)||_1')
%         xlabel('Iterations')
%     end

    % z-update
    [temp, ~] = radon3d_forward(LLL, fft(real(s),nfft,1), nt, ilow, ihigh);
    [~, temp] = radon3d_forward_adjoint(LLL, temp-D, nt, ilow, ihigh);
    z = s - step_size * temp;
    % m-update
    m_prev = m;
    m = wthresh(z,'s',step_size * reg_param);
    % q-update
    q_new = 0.5 * (1 + sqrt(1+4*(q_t^2)));
    % s-update
    s = m + (q_t-1)/q_new * (m - m_prev);
    q_t = q_new;

    pause(.0002)

end










