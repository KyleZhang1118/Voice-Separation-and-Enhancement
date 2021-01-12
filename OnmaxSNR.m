function [Y,SetupStruc] = OnmaxSNR(S,Transfer,SetupStruc)
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
covMa = cal_covMa4Block(SetupStruc.Block,blockWin);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
theta = 10^-3;
for i = 2:K_m
    X_f = permute(X(i,:,:),[2 3 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%maxSNR
    W_f = zeros(N,Num);
    for j = 1:Num
        R_d = covMa(:,:,i,j);
        R_i = zeros(N);
        for k = 1:Num
            if(k==j) continue;end
            R_i = R_i+covMa(:,:,i,k);
        end
        if(rcond(R_i)<theta)
            R_i = R_i+eye(N)*min(diag(R_i))*10^-3;
        end
        [E,~] = PCA(R_i\R_d,1,1);
        W_f(:,j) = E*sqrt(E'*R_i*R_i*E/N)/(E'*R_i*E);  %%%%%% WAN post filtering
    end
    W_f = W_f';
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