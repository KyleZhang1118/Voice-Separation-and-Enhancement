function [Y,W,SetupStruc] = Process_FastMNMF2(s,Transfer,SetupStruc)
K = SetupStruc.FastMNMF2.K;
hop = SetupStruc.FastMNMF2.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.FastMNMF2.win = win;  % Preserve 'win' in 'SetupStruc'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(s,2);
% s = s/max(max(s));
for i = 1:N
    X(:,:,i) = fft(enframe(s(:,i),win,hop)');
end
frame_N = size(X,2);
K_m = K/2+1;
Num = size(Transfer,3);
Y = zeros((frame_N-1)*hop+K,Num);
Y_f = zeros(Num,frame_N,K);
%%%%%%%%%%%%%%%%%%%%%%%%%% First initialization
epsi = 1e-7;
L = 2;        %%%%% the initial number of NMF basis
X_sp = zeros(K_m,frame_N,N);
Y_sp = zeros(K_m,frame_N,N);
T = max(rand(K_m,L,Num),epsi);
V = max(rand(L,frame_N,Num),epsi);
G = max(eye(Num,N),0.01);
if(N>Num)
    i = Num+1;
    j = 1;
    while(i<=N)
        G(j,i) = 1;
        i = i+1;
        j = j+1;
        if(j>Num)
            j = 1;
        end
    end
end
G = G./repmat(sum(G,2),[1,N]);
G = G/max(sum(G));
Q = eye(N);
Q = repmat(Q,1,1,K_m);
theta = 10^-6;
X_Norm = X;
for i = 1:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialization
    X_f = X_f./repmat(max(max(abs(X_f),[],2),epsi),[1,frame_N]);
    X_Norm(i,:,:) = X_f.';
    X_temp = abs(Q(:,:,i)*X_f).^2;
    X_sp(i,:,:) = X_temp'; 
end
for i = 1:N
    Y_temp = zeros(K_m,frame_N);
    for i_Num = 1:Num
        Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i);
    end
    Y_sp(:,:,i) = Y_temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Graudal iterations
% pObj = inf;
% A = zeros(1001,2)-1; %%%% Show the decrease of the value of cost funtion, ILRMA max iterations 1000
for iteration = 1:50
    %%%%% MU of NMF 
    for i = 1:Num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update T
        Gm = permute(repmat(G(i,:),K_m,1,frame_N),[1 3 2]);
        T(:,:,i) = max(T(:,:,i).*sqrt((sum(Gm.*X_sp.*Y_sp.^(-2),3)*V(:,:,i)')./max(sum(Gm./Y_sp,3)*V(:,:,i)',epsi)),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i_N);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update V
        V(:,:,i) = max(V(:,:,i).*sqrt((T(:,:,i)'*sum(Gm.*X_sp.*Y_sp.^(-2),3))./max(T(:,:,i)'*sum(Gm./Y_sp,3),epsi)),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i_N);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update G
        G(i,:) = max(G(i,:).*permute(sqrt(sum(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N).*X_sp.*Y_sp.^(-2)))./max(sum(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N)./Y_sp)),epsi)),[1 3 2]),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i_N);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
    end
    %%%%% IP of AuxIVA
%     dlw = 0;
    for i = 1:K_m
        X_f = permute(X_Norm(i,:,:),[3 2 1]);
%         dlw = dlw +log(abs(det(Q(:,:,i)))+epsi);
        %%%%% AuxIVA
        for i_n = 1:N
            G_ = permute(Y_sp(i,:,i_n),[3 2 1]);
            G_ = repmat(G_,N,1);
            Vk = (X_f./(G_+epsi))*X_f'/frame_N;
            if rcond(Vk)<theta
                Vk = Vk+eye(N)*max(eig(Vk))*theta;
            end
            wk = inv(Q(:,:,i)*Vk);
            wk = wk(:,i_n);
            wk = wk/(sqrt(wk'*Vk*wk)+epsi);
            Q(i_n,:,i) = wk';
        end
        X_temp = abs(Q(:,:,i)*X_f).^2;
        X_sp(i,:,:) = X_temp'; 
    end
%     Obj = ((sum(sum(sum(log(Y_sp+epsi))))+sum(sum(sum(X_sp./(Y_sp+epsi)))))/frame_N-2*dlw)/(Num*K_m);
%     dObj = pObj-Obj;
%     pObj = Obj;
%     A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
%     if(abs(dObj)/abs(Obj)<theta)
%        break;
%     end
    %%%%%%%% Adjust the scales
    for i = 1:K_m
        miu = trace(Q(:,:,i)*Q(:,:,i)')/N;
        Q(:,:,i) = Q(:,:,i)/sqrt(miu);
        T(i,:,:) = T(i,:,:)/miu;
    end
    for i = 1:Num
        phi = sum(G(i,:));
        G(i,:) = G(i,:)/phi;
        T(:,:,i) = T(:,:,i)*phi;
    end
    vv = sum(T,1);
    T = T./repmat(vv,K_m,1,1);
    V = V.*repmat(permute(vv,[2 1 3]),1,frame_N,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Iteration with larger L
L = 16;        %%%%% the initial number of NMF basis
Y_sp = zeros(K_m,frame_N,N);
T = max(rand(K_m,L,Num),epsi);
V = max(rand(L,frame_N,Num),epsi);
for i = 1:N
    Y_temp = zeros(K_m,frame_N);
    for i_Num = 1:Num
        Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i);
    end
    Y_sp(:,:,i) = Y_temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Graudal iterations
max_iteration = 150;
pObj = inf;
A = zeros(1001,2)-1; %%%% Show the decrease of the value of cost funtion, ILRMA max iterations 1000
for iteration = 1:max_iteration
    %%%%% MU of NMF 
    for i = 1:Num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update T
        Gm = permute(repmat(G(i,:),K_m,1,frame_N),[1 3 2]);
        T(:,:,i) = max(T(:,:,i).*sqrt((sum(Gm.*X_sp.*Y_sp.^(-2),3)*V(:,:,i)')./max(sum(Gm./Y_sp,3)*V(:,:,i)',epsi)),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i_N);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update V
        V(:,:,i) = max(V(:,:,i).*sqrt((T(:,:,i)'*sum(Gm.*X_sp.*Y_sp.^(-2),3))./max(T(:,:,i)'*sum(Gm./Y_sp,3),epsi)),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i_N);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update G
        G(i,:) = max(G(i,:).*permute(sqrt(sum(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N).*X_sp.*Y_sp.^(-2)))./max(sum(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N)./Y_sp)),epsi)),[1 3 2]),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    Y_temp = Y_temp+T(:,:,i_Num)*V(:,:,i_Num)*G(i_Num,i_N);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
    end
    %%%%% IP of AuxIVA
    dlw = 0;
    for i = 1:K_m
        X_f = permute(X_Norm(i,:,:),[3 2 1]);
        dlw = dlw +log(abs(det(Q(:,:,i)))+epsi);
        %%%%% AuxIVA
        for i_n = 1:N
            G_ = permute(Y_sp(i,:,i_n),[3 2 1]);
            G_ = repmat(G_,N,1);
            Vk = (X_f./(G_+epsi))*X_f'/frame_N;
            if rcond(Vk)<theta
                Vk = Vk+eye(N)*max(eig(Vk))*theta;
            end
            wk = inv(Q(:,:,i)*Vk);
            wk = wk(:,i_n);
            wk = wk/(sqrt(wk'*Vk*wk)+epsi);
            Q(i_n,:,i) = wk';
        end
        X_temp = abs(Q(:,:,i)*X_f).^2;
        X_sp(i,:,:) = X_temp'; 
    end
    Obj = ((sum(sum(sum(log(Y_sp+epsi))))+sum(sum(sum(X_sp./(Y_sp+epsi)))))/frame_N-2*dlw)/(Num*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
%     if(abs(dObj)/abs(Obj)<theta)
%        break;
%     end
    %%%%%%% Adjust the scales
    for i = 1:K_m
        miu = trace(Q(:,:,i)*Q(:,:,i)')/N;
        Q(:,:,i) = Q(:,:,i)/sqrt(miu);
        T(i,:,:) = T(i,:,:)/miu;
    end
    for i = 1:Num
        phi = sum(G(i,:));
        G(i,:) = G(i,:)/phi;
        T(:,:,i) = T(:,:,i)*phi;
    end
    vv = sum(T,1);
    T = T./repmat(vv,K_m,1,1);
    V = V.*repmat(permute(vv,[2 1 3]),1,frame_N,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Post processing
% W = conj(permute(Transfer,[3 2 1]))/N;
W = 0;
lamda = zeros(K_m,frame_N,Num);
for i = 1:Num
    lamda(:,:,i) = T(:,:,i)*V(:,:,i);
end
lamdag = zeros(Num,N);
for i = 1:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    Q_i = Q(:,:,i);
    if rcond(Q_i)<theta
        Q_i = Q_i+eye(N)*max(eig(Q_i))*theta;
    end
    Q_ii = inv(Q_i);
    for i_T = 1:frame_N
        for i_N = 1:Num
            lamdag(i_N,:) = lamda(i,i_T,i_N)*G(i_N,:);
        end
        for i_N = 1:Num
            Y_f(i_N,i_T,i) = Q_ii(1,:)*diag(lamdag(i_N,:)./sum(lamdag))*Q(:,:,i)*X_f(:,i_T);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(i~=K_m && i~=1)
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