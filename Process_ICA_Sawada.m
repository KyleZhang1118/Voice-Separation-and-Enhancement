function [Y,W,SetupStruc] = Process_ICA_Sawada(s,Transfer,SetupStruc)
K = SetupStruc.ICA_Sawada.K;
hop = SetupStruc.ICA_Sawada.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.ICA_Sawada.win = win;  % Preserve 'win' in 'SetupStruc'
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
Y_P = zeros(frame_N,Num,K_m);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
theta = 10^-4;
W = zeros(Num,N,K_m);
A = zeros(1001,K/2)-1; %%%% Show the decrease of non-linear correlation, ICA max iterations 1000
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PCA and ICA processing
    [E,D] = PCA(X_f,1,Num);
    V = sqrt(D)\E';
    X_f = V*X_f;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Steer = permute(Transfer(i,:,:),[2 3 1]);
%     Ori = V*Steer;
%     if rcond(Ori)<theta
%         Ori = Ori+eye(Num)*min(diag(Ori));
%     end
%     [Y_,W_ICA,A] = FDICA(X_f,inv(Ori),A,i);  %%% 'A', 'i' record the decrease for observation    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Y_,W_ICA,A] = FDICA(X_f,eye(Num),A,i);  %%% 'A', 'i' record the decrease for observation      
    W(:,:,i) = W_ICA*V;
    Y_P(:,:,i) = Y_.';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process the ambiguity of permutation and amplitude
P = Permu_Sawada(W,Y_P,SetupStruc,'all');  %%%% Options: 'DOA','cor', 'all'
for i = 2:K_m
    W(:,:,i) = P(:,:,i)*W(:,:,i);
    Y_ = permute(Y_P(:,:,i),[2 1 3]);
    Y_ = P(:,:,i)*Y_;
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