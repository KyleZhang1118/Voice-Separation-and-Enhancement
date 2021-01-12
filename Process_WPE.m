function [Y,SetupStruc] = Process_WPE(s,SetupStruc)
K = SetupStruc.WPE.K;
hop = SetupStruc.WPE.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.WPE.win = win;  % Preserve 'win' in 'SetupStruc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(s,2);
for i = 1:N
    X(:,:,i) = fft(enframe(s(:,i),win,hop)');
end
frame_N = size(X,2);
K_m = K/2+1;
Y = zeros((frame_N-1)*hop+K,N);
Y_f = zeros(size(X));
%%%%%%%%%%%%%%%%%%%%%%%%%%WPEÈ¥»ìÏì
Y_f(1:K_m,:,:) = WPE_DM(X(1:K_m,:,:));  %%%% 
Y_f(K_m+1:end,:,:) = conj(Y_f(K_m-1:-1:2,:,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:N
    Y(:,i) = overlapadd(real(ifft(Y_f(:,:,i)))',win,hop);
end

return;