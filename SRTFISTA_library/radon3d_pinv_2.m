function [M, m] = radon3d_pinv_2(LLL, LLLINV, DD1, DD2, nt, ilow, ihigh)
%forward radon in frequency domain -- return result in frequency and time
%domain @Dr. Olugboji

[nfft, ~, nq] = size(LLL);
M = zeros(nfft, nq);


for ifreq=ilow:ihigh
    
    L = squeeze( LLL(ifreq, :,:) );
    y = DD1(ifreq,:)';
    pinv_L = squeeze( LLLINV(ifreq, :,:) );
    
    xa = L'*y + (DD2(ifreq, :)');

    x  =  pinv_L \ xa;   % do pseudo inverse here
    
    M(ifreq,:) = x';
    M(nfft+2-ifreq,:) = conj(x)';
end

 M(nfft/2+1,:) = zeros(1,nq);
 m = real(ifft(M,[],1));
 m = m(1:nt,:);

return