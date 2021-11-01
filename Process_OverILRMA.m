function [Y,W,SetupStruc] = Process_OverILRMA(s,Transfer,SetupStruc)
K = SetupStruc.OverILRMA.K;
hop = SetupStruc.OverILRMA.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.OverILRMA.win = win;  % Preserve 'win' in 'SetupStruc'
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
X_sp = zeros(N,frame_N,K_m);
W_IVA = zeros(Num,N,K_m);
W_Over = zeros(N,N,K_m);
V_sp = zeros(N,N,K_m);
theta = 10^-6;
X = permute(X,[3 2 1]);
E1 = [eye(Num) zeros(Num,N-Num)];
E2 = [zeros(N-Num,Num) eye(N-Num)];
for i = 1:K_m
    X_f = X(:,:,i);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialization
%     W_o = [eye(Num) zeros(Num,N-Num)];
%     [E,D] = PCA(X_f,1,Num);
    V = diag(max(abs(X_f),[],2).^(-1));
    W = eye(Num,N);
    y_f = W*V*X_f;
    V_sp(:,:,i) = V;
    X_sp(:,:,i) = V*X_f;
    W_IVA(:,:,i) = W;
    Cf = X_f*X_f'/frame_N;
    Jf = E1*Cf*W';
    if rcond(Jf)<theta
       Jf = Jf+eye(Num)*max(eig(Jf))*theta;
    end
    Jf = E2*Cf*W'/Jf;
    W_ = [W;[Jf -eye(N-Num)]];
    W_Over(:,:,i) = W_;
%     W_Over(:,:,i) = [W_o;[zeros(N-Num,Num) -eye(N-Num)]];
    Y_f(:,:,i) = y_f;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% IVA iterations
max_iteration = 200;
epsi = 1e-7;
L = 2;  %%%%% the number of NMF basis
T = max(rand(K_m,L,Num),epsi);
V = max(rand(L,frame_N,Num),epsi);
R = zeros(K_m,frame_N,Num);
for i = 1:Num
    R(:,:,i) = T(:,:,i)*V(:,:,i);
end
P = permute(Y_f(:,:,1:K_m),[3,2,1]);
P = abs(P).^2;
pObj = inf;
A = zeros(1001,2)-1; %%%% Show the decrease of non-linear correlation, IVA max iterations 1000
for iteration = 1:max_iteration
    %%%%% NMF
    for i = 1:Num
        T(:,:,i) = max(T(:,:,i).*sqrt(P(:,:,i).*R(:,:,i).^(-2)*V(:,:,i)'./(R(:,:,i).^(-1)*V(:,:,i)')),epsi);
        R(:,:,i) = T(:,:,i)*V(:,:,i);
        V(:,:,i) = max(V(:,:,i).*sqrt(T(:,:,i)'*(P(:,:,i).*R(:,:,i).^(-2))./(T(:,:,i)'*R(:,:,i).^(-1))),epsi);
        R(:,:,i) = T(:,:,i)*V(:,:,i);
    end
    dlw = 0;
    for i = 1:K_m
        W = W_IVA(:,:,i);
        W_ = W_Over(:,:,i);
        X_f = X_sp(:,:,i);
        Cf = X_f*X_f'/frame_N;
        dlw = dlw +log(abs(det(W_))+epsi);
        for i_n = 1:Num
            G_ = permute(R(i,:,i_n),[3 2 1]);
            G_ = repmat(G_+epsi,N,1);
            Vk = (X_f./G_)*X_f'/frame_N;
            if rcond(Vk)<theta
%                 coe = sort(eig(Vk),'descend');
                Vk = Vk+eye(N)*max(eig(Vk))*theta;
            end
            wk = inv(W_*Vk);
            wk = wk(:,i_n);
            wk = wk/(sqrt(wk'*Vk*wk)+epsi);
            W(i_n,:) = wk';
            Jf = E1*Cf*W';
            if rcond(Jf)<theta
                Jf = Jf+eye(Num)*max(eig(Jf))*theta;
            end
            Jf = E2*Cf*W'/Jf;
            W_ = [W;[Jf -eye(N-Num)]];
        end
        W_IVA(:,:,i) = W;
        W_Over(:,:,i) = W_;
        Y_f(:,:,i) = W*X_f;
    end
    P = permute(Y_f(:,:,1:K_m),[3,2,1]);
    P = abs(P).^2;
%     Obj = (sum(sum(Y_k))/frame_N-dlw)/(Num*K_m);
    Obj = ((sum(sum(sum(log(R+epsi))))+sum(sum(sum(P./R))))/frame_N-2*dlw)/(Num*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
    if(abs(dObj)/abs(Obj)<epsi)
        break;
    end
    for i = 1:Num
        lamda = sqrt(sum(sum(P(:,:,i))/(frame_N*K_m)));
        W_IVA(i,:,:) = W_IVA(i,:,:)/lamda;
        P(:,:,i) = P(:,:,i)/lamda^2;
        R(:,:,i) = R(:,:,i)/lamda^2;
        T(:,:,i) = T(:,:,i)/lamda^2;
    end
    W_Over(1:Num,:,:) = W_IVA;      
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