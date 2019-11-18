function [Y,SetupStruc] = OnDSB_Mask(S,Transfer,SetupStruc)
K = SetupStruc.K;
hop = SetupStruc.hop;
win = SetupStruc.win;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(S,2);
frame_N = size(S,3);
blockWin = repmat(win,[1 N frame_N]);
X = fft(S.*blockWin);
K_m = K/2+1;
Num = size(Transfer,3);
Y = zeros((frame_N-1)*hop+K,Num);
Y_f = zeros(K,frame_N,Num);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
W = conj(permute(Transfer,[3 2 1]))/N;
for i = 2:K_m
    X_f = permute(X(i,:,:),[2 3 1]);
    W_f = permute(W(:,:,i),[1 2 3]);
    Y_ = W_f*X_f;
    [Y_max,I] = max(Y_);
    Y_ = zeros(size(Y_));
    Y_(I) = Y_max;
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