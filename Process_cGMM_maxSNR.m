function [Y,W,SetupStruc] = Process_cGMM_maxSNR(s,Transfer,SetupStruc)
K = SetupStruc.cGMM_maxSNR.K;
hop = SetupStruc.cGMM_maxSNR.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.cGMM_maxSNR.win = win;  % Preserve 'win' in 'SetupStruc'
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
%covMa = cal_covMa(s,K,hop,win);
R_X = zeros(N,N,frame_N);
if(~isfield(SetupStruc.cGMM_maxSNR,'lamda'))
    R_initial = cal_covMa(SetupStruc.unS,K,hop,win);
    [lamda,R,SetupStruc.cGMM_maxSNR.Q] = cGMM(X(1:K_m,:,:),Transfer,20,R_initial);
    SetupStruc.cGMM_maxSNR.lamda = lamda;
    SetupStruc.cGMM_maxSNR.R = R;
else
    lamda = SetupStruc.cGMM_maxSNR.lamda;
    R = SetupStruc.cGMM_maxSNR.R;
end
Num_v = size(lamda,2);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
theta = 10^-4;
W = zeros(Num,N,K_m);
for i = 2:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    covMa = zeros(N,N,Num_v);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%way1
    for t = 1:frame_N
        R_X(:,:,t) = X_f(:,t)*X_f(:,t)';
    end
    for j = 1:Num_v
        R_temp = zeros(N);
        for t = 1:frame_N
            R_temp = R_temp+lamda(t,j,i)*R_X(:,:,t);
        end
        covMa(:,:,j) = R_temp/frame_N;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%way2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%maxSNR
    W_f = zeros(N,Num);
    for j = 1:Num
        R_d = covMa(:,:,j);
        R_i = zeros(N);
        for k = 1:Num_v
            if(k==j) continue;end
            R_i = R_i+covMa(:,:,k);
        end
        if(rcond(R_i)<theta)
            R_i = R_i+eye(N)*max(diag(R_i))*theta;
        end
        [E,~] = PCA(R_i\R_d,1,1);
        W_f(:,j) = E*sqrt(E'*R_i*R_i*E/N)/(E'*R_i*E);  %%%%%% WAN post filtering
    end
    W_f = W_f';
    W(:,:,i) = W_f;
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