function [Y,SetupStruc] = OnLCMV(S,Transfer,SetupStruc)
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
%%%%%%%%%%%%%%%%%%%%%%%%%% Initial the coefficients stored in 'SetupStruc'
if(~isfield(SetupStruc,'W'))
    SetupStruc.W = permute(Transfer,[2 3 1])/N;
    SetupStruc.V = zeros(N,Num,K_m);
    sign_Int = 1;
    SetupStruc.R = zeros(N,N,K_m);
    SetupStruc.WNG = [];
%     SetupStruc.miu = [];
%     SetupStruc.mag = [];
else
    sign_Int = 0;
end
WNG_set = 3;
WNG = zeros(1,K_m);
% miu_ = zeros(1,K_m);
% mag = zeros(1,K_m);
delta2 = 10^(WNG_set/10);
for i = 2:K_m
    X_f = permute(X(i,:,:),[2 3 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MVDR processing
    W_f = zeros(N,Num);
    Steer = permute(Transfer(i,:,:),[2 3 1]);
    if(sign_Int==1 || SetupStruc.alpha==1)
        R = X_f*X_f'/frame_N;       
    else
        R = (1-SetupStruc.alpha)*SetupStruc.R(:,:,i)+SetupStruc.alpha*(X_f*X_f')/frame_N;
        SetupStruc.R(:,:,i) = R;
    end
    miu = 2/3/sum(real(diag(R)));
    Wc0 = Steer/(Steer'*Steer);
    Pc0 = eye(N)-Steer/(Steer'*Steer)*Steer';
    for j = 1:Num
        Wc = Wc0(:,j);
        V = SetupStruc.V(:,j,i);
        W = SetupStruc.W(:,j,i);
        V = Pc0*(V-miu*R*W);  
        b2 = 1/delta2-1/(Wc'*Wc);
        if(V'*V<=b2)
            W = Wc+V;
        else
            W = Wc+sqrt(b2)*V/sqrt(V'*V);
        end
        W_f(:,j) = W;
        SetupStruc.V(:,j,i) = V;
    end
    SetupStruc.W(:,:,i) = W_f;
    WNG(i) = -10*log10(W_f(:,1)'*W_f(:,1));
%     mag(i) = abs(W_f(:,1)'*Steer(:,1));
%     miu_(i) = miu;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y_ = W_f'*X_f;
    Y_f(i,:,:) = Y_.';
    if(i~=K_m)
        Y_f(K+2-i,:,:) = Y_';
    end
end
SetupStruc.WNG = [SetupStruc.WNG;WNG];
% SetupStruc.miu = [SetupStruc.miu;miu_];
% SetupStruc.mag = [SetupStruc.mag;mag];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num
    Y(:,i) = overlapadd(real(ifft(Y_f(:,:,i)))',win,hop);
end
return;