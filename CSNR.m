function SNR = CSNR(S,order)
%%%% the 1st column is the signal, the others are the noises
N = S;
N(:,order) = zeros(size(N,1),1);
N = sum(N,2);
SNR = 10*log10(sum(S(:,order).^2)/sum(N.^2));
return