function [Y,W,SetupStruc] = Process_ICA_initial(s,Transfer,SetupStruc)
K = SetupStruc.ICA_initial.K;
hop = SetupStruc.ICA_initial.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.ICA_initial.win = win;  % Preserve 'win' in 'SetupStruc'
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
theta = 10^-4;
cond = zeros(K_m-1,2);
W = zeros(Num,N,K_m);
A = zeros(1001,K/2)-1; %%%% Show the decrease of non-linear correlation, ICA max iterations 1000
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PCA and ICA processing
    Steer = permute(Transfer(i,:,:),[2 3 1]);
    [E,D] = PCA(X_f,1,Num);
    V = sqrt(D)\E';
    X_f = V*X_f;
    Ori = V*Steer;
    cond(i-1,1) = rcond(Ori);  %%%%%%%
    if rcond(Ori)<theta
        Ori = Ori+eye(Num)*min(diag(Ori));
    end
    [Y_,W_ICA,A] = FDICA(X_f,inv(Ori),A,i);  %%% 'A', 'i' record the decrease for observation   
    W_inv = pinv(W_ICA*V);
    for ii = 1:Num
        Y_(ii,:) = Y_(ii,:)*W_inv(1,ii);
        W_ICA(ii,:) = W_ICA(ii,:)*W_inv(1,ii);
    end
    W(:,:,i) = W_ICA*V;     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y_f(i,:,:) = Y_.';
    if(i~=K_m)
        Y_f(K+2-i,:,:) = Y_';
    end
    cond(i-1,2) = rcond(Ori); %%%%%%%
end
SetupStruc.ICA_initial.A = A;
% figure
% plot(cond(:,1))
% hold on 
% plot(cond(:,2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num
    Y(:,i) = overlapadd(real(ifft(Y_f(:,:,i)))',win,hop);
end
return;