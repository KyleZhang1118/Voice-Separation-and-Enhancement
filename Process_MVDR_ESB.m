function [Y,W,SetupStruc] = Process_MVDR_ESB(s,Transfer,SetupStruc)
K = SetupStruc.MVDR_ESB.K;
hop = SetupStruc.MVDR_ESB.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.MVDR_ESB.win = win;  % Preserve 'win' in 'SetupStruc'
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
WNG =zeros(K_m,Num);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
theta = 10^-4;
W = zeros(Num,N,K_m);
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MVDR processing
    W_f = zeros(N,Num);
    Steer = permute(Transfer(i,:,:),[2 3 1]);
    [E,D] = PCA(X_f,1,N);
    E = E(:,1:Num);
    if(Num<N)
        D = D(1:Num,1:Num)-mean(diag(D(Num+1:end,Num+1:end)))*eye(Num);
    end
    if(~isempty(find(diag(D)<theta, 1)))
        D = D+min(diag(D))*eye(Num);
%         D = D+theta*eye(Num);
    end
    R_inv = E/D*E';
    for j = 1:Num
        h = Steer(:,j);
        W_f(:,j) = R_inv*h/(h'*R_inv*h);
    end
    W_f = W_f';
    W(:,:,i) = W_f;   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WNG(i,:) = -10*log10(diag(W_f*W_f')');
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
% autoPlot(WNG,'MVDR_ESB');
return;