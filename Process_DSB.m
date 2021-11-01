function [Y,W,SetupStruc] = Process_DSB(s,Transfer,SetupStruc)
K = SetupStruc.DSB.K;
hop = SetupStruc.DSB.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.DSB.win = win;  % Preserve 'win' in 'SetupStruc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(s,2);
for i = 1:N
    X(:,:,i) = fft(enframe(s(:,i),win,hop)');
end
frame_N = size(X,2);
K_m = K/2+1;
Num = size(Transfer,3);
Y = zeros((frame_N-1)*hop+K,Num);
Y_f = zeros(size(X,1),size(X,2),Num);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
W = conj(permute(Transfer,[3 2 1]))/N;
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    W_f = permute(W(:,:,i),[1 2 3]);
    Y_ = W_f*X_f;
    Y_f(i,:,:) = Y_.';
    if(i~=K_m)
        Y_f(K+2-i,:,:) = Y_';
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num
    Y(:,i) = overlapadd(real(ifft(Y_f(:,:,i)))',win,hop);
end

return;