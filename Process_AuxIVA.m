function [Y,W,SetupStruc] = Process_AuxIVA(s,Transfer,SetupStruc)
K = SetupStruc.AuxIVA.K;
hop = SetupStruc.AuxIVA.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.AuxIVA.win = win;  % Preserve 'win' in 'SetupStruc'
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
X_sp = zeros(Num,frame_N,K_m);
W_IVA = zeros(Num,Num,K_m);
V_sp = zeros(Num,N,K_m);
theta = 10^-6;
for i = 1:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize W by PCA
    [E,D] = PCA(X_f,1,Num);
    V = sqrt(D)\E';
    V_sp(:,:,i) = V;
    X_sp(:,:,i) = V*X_f;
    %%%%%%%%%%%% Adjust amplitude of 'w'
    W_o = eye(Num);
    y_f = W_o*V*X_f;
%     norm = max(abs(y_f),[],2);
%     if(norm>10)
%         norm = repmat(norm,1,Num);
%         W_o = W_o./norm;
%         y_f = W_o*V*X_f;
%     end
    W_IVA(:,:,i) = W_o;
    Y_f(:,:,i) = y_f;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% IVA iterations
max_iteration = 200;
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
        X_f = X_sp(:,:,i);
        dlw = dlw +log(abs(det(W))+epsi);
        for i_n = 1:Num
            G_ = Y_k(i_n,:).^-1;
            G_ = repmat(G_,Num,1);
            Vk = (G_.*X_f)*X_f'/frame_N;
            if rcond(Vk)<theta
                Vk = Vk+eye(Num)*max(eig(Vk))*theta;
            end
            wk = inv(W*Vk);
            wk = wk(:,i_n);
            wk = wk/(sqrt(wk'*Vk*wk)+epsi);
            W(i_n,:) = wk';
        end
        W_IVA(:,:,i) = W;
        Y_f(:,:,i) = W*X_f;
    end
    Obj = (sum(sum(Y_k))/frame_N-2*dlw)/(Num*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
    if(abs(dObj)/abs(Obj)<theta)
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Post processing
W = zeros(Num,N,K_m);
Y_f(:,:,1) = zeros(Num,frame_N);
for i = 2:K_m
    W_inv = pinv(W_IVA(:,:,i)*V_sp(:,:,i));
    for ii = 1:Num
        Y_f(ii,:,i) = Y_f(ii,:,i)*W_inv(1,ii);
        W_IVA(ii,:,i) = W_IVA(ii,:,i)*W_inv(1,ii);
    end
    W(:,:,i) = W_IVA(:,:,i)*V_sp(:,:,i);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i~=K_m)
        Y_f(:,:,K+2-i) = conj(Y_f(:,:,i));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recover signals
if(K/hop==2)
    win = ones(K,1);
end
for i = 1:Num
    y_temp = permute(Y_f(i,:,:),[3 2 1]);
    Y(:,i) = overlapadd(real(ifft(y_temp))',win,hop);
end
return;