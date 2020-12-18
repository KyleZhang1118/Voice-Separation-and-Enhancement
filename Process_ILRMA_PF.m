function [Y,W,SetupStruc] = Process_ILRMA_PF(s,Transfer,SetupStruc)
K = SetupStruc.ILRMA_PF.K;
hop = SetupStruc.ILRMA_PF.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.ILRMA_PF.win = win;  % Preserve 'win' in 'SetupStruc'
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
W_ILRMA = zeros(Num,Num,K_m);
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
    W_ILRMA(:,:,i) = W_o;
    Y_f(:,:,i) = y_f;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ILRMA iterations
epsi = 1e-6;
max_iteration = 300;
L = 2*Num;  %%%%% the number of NMF basis
Z = max(rand(Num,L),epsi);
Z = Z./(ones(Num)*Z);
T = max(rand(K_m,L),epsi);
V = max(rand(L,frame_N),epsi);
R = zeros(K_m,frame_N,Num);
for i = 1:Num
    R(:,:,i) = (ones(K_m,1)*Z(i,:)).*T*V;
end
P = permute(Y_f(:,:,1:K_m),[3,2,1]);
P = abs(P).^2;
% K_sp = sqrt(K_m);
pObj = inf;
bZ = zeros(size(Z));
bT = zeros(size(T));
bV = zeros(size(V));
lamd = zeros(Num,1);
A = zeros(1001,2)-1; %%%% Show the decrease of the value of cost funtion, ILRMA max iterations 1000
for iteration = 1:max_iteration
    %%%%%% NMF with partioning function Z
    for i = 1:Num
        bZ(i,:) = sqrt((T'*(P(:,:,i).*(R(:,:,i).^(-2)))).*V*ones(frame_N,1)./((T'*R(:,:,i).^(-1)).*V*ones(frame_N,1)))';
    end
    Z = max(Z./(ones(Num)*Z),epsi);
    for i = 1:Num
        R(:,:,i) = (ones(K_m,1)*Z(i,:)).*T*V;
    end
    for i = 1:K_m
        P_temp = permute(P(i,:,:),[2 3 1]);
        R_temp = permute(R(i,:,:),[2 3 1]);
        bT(i,:) = sqrt((V*(P_temp.*(R_temp.^(-2)))).*Z'*ones(Num,1)./((V*R_temp.^(-1)).*Z'*ones(Num,1)))';
    end
    T = max(T.*bT,epsi);
    for i = 1:Num
        R(:,:,i) = (ones(K_m,1)*Z(i,:)).*T*V;
    end
    for i = 1:frame_N
        P_temp = permute(P(:,i,:),[1 3 2]);
        R_temp = permute(R(:,i,:),[1 3 2]);
        bV(:,i) = sqrt((T'*(P_temp.*(R_temp.^(-2)))).*Z'*ones(Num,1)./((T'*R_temp.^(-1)).*Z'*ones(Num,1)));
    end
    V = max(V.*bV,epsi);
    for i = 1:Num
        R(:,:,i) = (ones(K_m,1)*Z(i,:)).*T*V;
    end
    %%%%% AuxIVA
    dlw = 0;
    for i = 1:K_m
        W = W_ILRMA(:,:,i);
        X_f = X_sp(:,:,i);
       dlw = dlw +log(abs(det(W))+epsi);
        for i_n = 1:Num
            G_ = permute(R(i,:,i_n),[3 2 1]);
            G_ = repmat(sqrt(G_),Num,1);
            Vk = (X_f./G_)*X_f'/frame_N;
            if rcond(Vk)<theta
                Vk = Vk+eye(Num)*min(eig(Vk))*theta;
            end
            wk = inv(W*Vk);
            wk = wk(:,i_n);
            wk = wk/(sqrt(wk'*Vk*wk)+epsi);
            W(i_n,:) = wk';
        end
        W_ILRMA(:,:,i) = W;
        Y_f(:,:,i) = W*X_f;
    end
    P = permute(Y_f(:,:,1:K_m),[3,2,1]);
    P = abs(P).^2;
    Obj = ((sum(sum(sum(log(R+epsi))))+sum(sum(sum(P./R))))/frame_N-2*dlw)/(Num*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
    if(abs(dObj)/abs(Obj)<theta)
       break;
    end
    for i = 1:Num
        lamda = sqrt(sum(sum(P(:,:,i))/(frame_N*K_m)));
        lamd(i) = lamda;
        for i_k = 1:K_m
            W_ILRMA(i,:,i_k) = W_ILRMA(i,:,i_k)/lamda;
            for l = 1:L
                T(i_k,l) = T(i_k,l)*sum(Z(:,l))/lamda^2;
            end
        end
        P(:,:,i) = P(:,:,i)/lamda^2;
        R(:,:,i) = R(:,:,i)/lamda^2;
    end
    for l = 1:L
        Z_temp = sum(Z(:,l)./lamd.^(2));
        for i = 1:Num
            Z(i,l) = Z(i,l)/lamd(i)^2/Z_temp;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Post processing
W = zeros(Num,N,K_m);
Y_f(:,:,1) = zeros(Num,frame_N);
for i = 2:K_m
    W_inv = pinv(W_ILRMA(:,:,i)*V_sp(:,:,i));
    for ii = 1:Num
        Y_f(ii,:,i) = Y_f(ii,:,i)*W_inv(1,ii);
        W_ILRMA(ii,:,i) = W_ILRMA(ii,:,i)*W_inv(1,ii);
    end
    W(:,:,i) = W_ILRMA(:,:,i)*V_sp(:,:,i);     
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