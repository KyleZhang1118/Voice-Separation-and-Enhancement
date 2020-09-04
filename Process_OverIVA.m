function [Y,W,SetupStruc] = Process_OverIVA(s,Transfer,SetupStruc)
K = SetupStruc.OverIVA.K;
hop = SetupStruc.OverIVA.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.OverIVA.win = win;  % Preserve 'win' in 'SetupStruc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(s,2);
for i = 1:N
    X(:,:,i) = fft(enframe(s(:,i),win,hop)');
end
frame_N = size(X,2);
K_m = K/2+1;
Num = size(Transfer,3);
Y = zeros((frame_N-1)*hop+K,Num);
Y_f = zeros(Num,frame_N,K);
%%%%%%%%%%%%%%%%%%%%%%%%%% Obtain processing matrix 'W'
W_IVA = zeros(Num,N,K_m);
W_Over = zeros(N,N,K_m);
theta = 10^-4;
X = permute(X,[3 2 1]);
E1 = [eye(Num) zeros(Num,N-Num)];
E2 = [zeros(N-Num,Num) eye(N-Num)];
for i = 1:K_m
    X_f = X(:,:,i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialization
%     W_o = [eye(Num) zeros(Num,N-Num)];
    [E,D] = PCA(X_f,1,Num);
    W = sqrt(D)\E';
    y_f = W*X_f;
    W_IVA(:,:,i) = W;
    Cf = X_f*X_f'/frame_N;
    Jf = E1*Cf*W';
    if rcond(Jf)<theta
       Jf = Jf+eye(Num)*min(eig(Jf))*theta;
    end
    Jf = E2*Cf*W'/Jf;
    W_ = [W;[Jf -eye(N-Num)]];
    W_Over(:,:,i) = W_;
%     W_Over(:,:,i) = [W_o;[zeros(N-Num,Num) -eye(N-Num)]];
    Y_f(:,:,i) = y_f;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% IVA iterations
max_iteration = 100;
Y_k = zeros(Num,frame_N);
epsi = 1e-6;
pObj = inf;
A = zeros(1001,2)-1; %%%% Show the decrease of non-linear correlation, IVA max iterations 1000
for iteration = 1:max_iteration
    for i = 1:Num
        y_temp = permute(Y_f(i,:,:),[3 2 1]);
        Y_k(i,:) = sqrt(sum(abs(y_temp(1:K_m,:)).^2))+epsi;
    end
    dlw = 0;
    for i = 1:K_m
        W = W_IVA(:,:,i);
        W_ = W_Over(:,:,i);
        X_f = X(:,:,i);
        Cf = X_f*X_f'/frame_N;
        dlw = dlw +log(abs(det(W_))+epsi);
        for i_n = 1:Num
            G_ = Y_k(i_n,:).^-1;
            G_ = repmat(G_,N,1);
            Vk = (G_.*X_f)*X_f'/frame_N;
            if rcond(Vk)<theta
                coe = sort(eig(Vk),'descend');
                Vk = Vk+eye(N)*coe(Num)*theta;
            end
            wk = inv(W_*Vk);
            wk = wk(:,i_n);
            wk = wk/(sqrt(wk'*Vk*wk)+epsi);
            W(i_n,:) = wk';
            Jf = E1*Cf*W';
            if rcond(Jf)<theta
                Jf = Jf+eye(Num)*min(eig(Jf))*theta;
            end
            Jf = E2*Cf*W'/Jf;
            W_ = [W;[Jf -eye(N-Num)]];
        end
        W_IVA(:,:,i) = W;
        Y_f(:,:,i) = W*X_f;
    end
    Obj = (sum(sum(Y_k))/frame_N-dlw)/(Num*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
%     if(abs(dObj)/abs(Obj)<epsi)
%         break;
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Post processing
% W = zeros(Num,N,K_m);
Y_f(:,:,1) = zeros(Num,frame_N);
for i = 2:K_m
    W_inv = pinv(W_IVA(:,:,i));
    for ii = 1:Num
        Y_f(ii,:,i) = Y_f(ii,:,i)*W_inv(1,ii);
        W_IVA(ii,:,i) = W_IVA(ii,:,i)*W_inv(1,ii);
    end
%     W(:,:,i) = W_IVA(:,:,i);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i~=K_m)
        Y_f(:,:,K+2-i) = conj(Y_f(:,:,i));
    end
end
W = W_IVA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num
    y_temp = permute(Y_f(i,:,:),[3 2 1]);
    Y(:,i) = overlapadd(real(ifft(y_temp))',win,hop);
end
return;