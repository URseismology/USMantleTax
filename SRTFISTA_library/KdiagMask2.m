function [Kmask, qps_q] = KdiagMask2(M, taus, qs, H, Vp, Vs, target, time_range, amplitude_option)
% M: input matrix
% taus: time array
% qs: q values array
% H, Vp, Vs, Velocity model
% time_range: [t_start, t_end], desired time window to keep
% amplitude_option: 'positive', 'negative' or 'both' to keep respective values

nn = 100;
hbot = linspace(20, 200, nn);
tps_umd = zeros(nn,1);
qps_umd = zeros(nn, 1);

for in = 1:nn
    H(end-1) = hbot(in);
    [q, tau] = get_q_t(H, Vp, Vs);

    tps_umd(in) = tau(target, 1);
    qps_umd(in) = q(target, 1);
end

qps_q = interp1(tps_umd, qps_umd, taus, "linear", "extrap");

[nq, nt] = size(M);

Kmask = ones(nq, nt);

for it = 1:nt

    qlq = qps_q(it) - 20; qlq(qlq<0) = 0;
    qhq = qps_q(it) + 600;

    ilwq = qs < (qlq);
    ihiq = qs > (qhq);

    Kmask(ilwq, it) = 0;
    Kmask(ihiq, it) = 0;

end

% Mask based on amplitude option
if strcmp(amplitude_option, 'positive')
    Kmask(M <= 0) = 0;
elseif strcmp(amplitude_option, 'negative')
    Kmask(M >= 0) = 0;
% 'both' option will not change the mask
end

% Mask everything outside defined time range
Kmask(:, taus < time_range(1) | taus > time_range(2)) = 0;

end
