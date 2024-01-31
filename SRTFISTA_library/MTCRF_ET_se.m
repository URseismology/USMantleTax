function [time, radRF, transRF, varR, varT] = ...
    MTCRF_ET_se(R, T, Z, dnoi, PSI, K, Dt, taper_win, olap)

N = length(Z); % length of input/output
Ntaper = round(taper_win / Dt); % length of each taper in samples
ofac   = 1 / (1 - olap); % overlap factor
Nwin   = round(ofac * N / Ntaper); % number of moving windows

YZ = repmat(complex(0), N, K);
YR = repmat(complex(0), N, K);
YT = repmat(complex(0), N, K);
So = repmat(complex(0), N, K);

for j = 1:Nwin

    js = round(j * Ntaper / ofac);
    je = round(js + Ntaper); % start and end indices for this window

    if je > N
        break;
    end

    for k = 1:K

        w = zeros(N, 1);
        w(js:je - 1) = PSI(:, k); % k-th slepian taper

        % taper the data and do FFT
        Dat = w .* [Z, R, T, dnoi];
        fDat = fft(Dat);

        % k-th multi-taper spectral estimates and cumulate
        YZ(:,k) = YZ(:,k) + fDat(:,1);
        YR(:,k) = YR(:,k) + fDat(:,2);
        YT(:,k) = YT(:,k) + fDat(:,3);
        So(:,k) = So(:,k) + fDat(:,4);

    end
end

% cross-spectrum between components

SZT = repmat(complex(0), N, 1);
SZR = repmat(complex(0), N, 1);
SZZ = repmat(complex(0), N, 1);
SRR = repmat(complex(0), N, 1);
STT = repmat(complex(0), N, 1);
So  = repmat(complex(0), N, 1);

for k = 1:K

    SZR = SZR + conj(YZ(:,k)) .* YR(:,k);
    SZT = SZT + conj(YZ(:,k)) .* YT(:,k);
    SZZ = SZZ + conj(YZ(:,k)) .* YZ(:,k);
    SRR = SRR + conj(YR(:,k)) .* YR(:,k);
    STT = STT + conj(YT(:,k)) .* YT(:,k);

end

% frequency domain RFs through deconvolution (devision)

HR = SZR ./ (SZZ + So);
HT = SZT ./ (SZZ + So);

% Coherence estimates
CR = abs(SZR) ./ (sqrt(SZZ) .* sqrt(SRR) );
CT = abs(SZT) ./ (sqrt(SZZ) .* sqrt(STT) );

% variance

varR = ((1 - CR.^2) ./ ( (K-1) .*  (CR.^2))) .* (abs(HR).^2) ;
varT = ((1 - CT.^2) ./ ( (K-1) .*  (CT.^2))) .* (abs(HT).^2) ;

% return frequency domain RFs

radRF = HR;
transRF = HT;
time = ((0:N-1)-floor(N/2)) * Dt;