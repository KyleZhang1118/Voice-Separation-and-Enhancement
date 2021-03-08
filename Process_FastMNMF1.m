function [Y,W,SetupStruc] = Process_FastMNMF1(s,Transfer,SetupStruc)
K = SetupStruc.FastMNMF1.K;
hop = SetupStruc.FastMNMF1.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.FastMNMF1.win = win;  % Preserve 'win' in 'SetupStruc'
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
%%%%%%%%%%%%%%%%%%%%%%%%%% First initialization
epsi = 1e-6;
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
G_scale = max(sum(G));
G = repmat(G,1,1,K_m);
Q = eye(N);
Q = repmat(Q,1,1,K_m);
theta = 10^-6;
X_Norm = X;
for i = 1:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialization
    X_f = X_f./repmat(max(max(abs(X_f),[],2),epsi),[1,frame_N]);
    X_f = X_f*G_scale;
    X_Norm(i,:,:) = X_f.';
    X_temp = abs(Q(:,:,i)*X_f).^2;
    X_sp(i,:,:) = X_temp'; 
end
for i = 1:N
    Y_temp = zeros(K_m,frame_N);
    for i_Num = 1:Num
        G_temp = permute(G(i_Num,i,:),[3 1 2]);
        TG = T(:,:,i_Num).*repmat(G_temp,1,L);
        Y_temp = Y_temp+TG*V(:,:,i_Num);
    end
    Y_sp(:,:,i) = Y_temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Graudal iterations
% max_iteration = 1000;
% pObj = inf;
% A = zeros(1001,2)-1; %%%% Show the decrease of the value of cost funtion, ILRMA max iterations 1000
for iteration = 1:50
    %%%%% MU of NMF 
    for i = 1:Num
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update T
        Gm = permute(G(i,:,:),[3 2 1]);
        Gm = repmat(Gm,1,1,frame_N);
        Gm = permute(Gm,[1 3 2]);
        T(:,:,i) = max(T(:,:,i).*sqrt((sum(Gm.*X_sp.*Y_sp.^(-2),3)*V(:,:,i)')./(sum(Gm./Y_sp,3)*V(:,:,i)')),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    G_temp = permute(G(i_Num,i_N,:),[3 1 2]);
                    TG = T(:,:,i_Num).*repmat(G_temp,1,L);
                    Y_temp = Y_temp+TG*V(:,:,i_Num);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update V
        V(:,:,i) = max(V(:,:,i).*sqrt((T(:,:,i)'*sum(Gm.*X_sp.*Y_sp.^(-2),3))./(T(:,:,i)'*sum(Gm./Y_sp,3))),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    G_temp = permute(G(i_Num,i_N,:),[3 1 2]);
                    TG = T(:,:,i_Num).*repmat(G_temp,1,L);
                    Y_temp = Y_temp+TG*V(:,:,i_Num);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update G
        Gn = permute(G(i,:,:),[2 3 1]);
        G(i,:,:) = max(Gn.*sqrt(permute(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N).*X_sp.*Y_sp.^(-2),2),[3 1 2])./permute(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N)./Y_sp,2),[3 1 2])),0.01);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    G_temp = permute(G(i_Num,i_N,:),[3 1 2]);
                    TG = T(:,:,i_Num).*repmat(G_temp,1,L);
                    Y_temp = Y_temp+TG*V(:,:,i_Num);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
    end
    %%%%% IP of AuxIVA
%     dlw = 0;
    for i = 1:K_m
        X_f = permute(X_Norm(i,:,:),[3 2 1]);
%         dlw = dlw +log(abs(det(Q(:,:,i)))+epsi);
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
%     Obj = ((sum(sum(sum(log(Y_sp+epsi))))+sum(sum(sum(X_sp./Y_sp))))/frame_N-2*dlw)/(Num*K_m);
%     dObj = pObj-Obj;
%     pObj = Obj;
%     A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
%     if(abs(dObj)/abs(Obj)<theta)
%        break;
%     end
    %%%%%%%% Adjust the scales
    for i = 1:K_m
        for i_N = 1:N
            miu = Q(i_N,:,i)*Q(i_N,:,i)';
            Q(i_N,:,i) = Q(i_N,:,i)/sqrt(miu);
            G(:,i_N,i) = G(:,i_N,i)/miu;
        end
        for i_N = 1:Num
            phi = sum(G(i_N,:,i));
            G(i_N,:,i) = G(i_N,:,i)/phi;
            T(i,:,i_N) = T(i,:,i_N)*phi;
        end
    end
    vv = sum(T,1);
    T = T./repmat(vv,K_m,1,1);
    V = V.*repmat(permute(vv,[2 1 3]),1,frame_N,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Iteration with larger L
L = 16;        %%%%% the initial number of NMF basis
Y_sp = zeros(K_m,frame_N,N);
T = max(rand(K_m,L,Num),epsi);
V = max(rand(L,frame_N,Num),epsi);
for i = 1:N
    Y_temp = zeros(K_m,frame_N);
    for i_Num = 1:Num
        G_temp = permute(G(i_Num,i,:),[3 1 2]);
        TG = T(:,:,i_Num).*repmat(G_temp,1,L);
        Y_temp = Y_temp+TG*V(:,:,i_Num);
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
        Gm = permute(G(i,:,:),[3 2 1]);
        Gm = repmat(Gm,1,1,frame_N);
        Gm = permute(Gm,[1 3 2]);
        T(:,:,i) = max(T(:,:,i).*sqrt((sum(Gm.*X_sp.*Y_sp.^(-2),3)*V(:,:,i)')./(sum(Gm./Y_sp,3)*V(:,:,i)')),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    G_temp = permute(G(i_Num,i_N,:),[3 1 2]);
                    TG = T(:,:,i_Num).*repmat(G_temp,1,L);
                    Y_temp = Y_temp+TG*V(:,:,i_Num);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update V
        V(:,:,i) = max(V(:,:,i).*sqrt((T(:,:,i)'*sum(Gm.*X_sp.*Y_sp.^(-2),3))./(T(:,:,i)'*sum(Gm./Y_sp,3))),epsi);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    G_temp = permute(G(i_Num,i_N,:),[3 1 2]);
                    TG = T(:,:,i_Num).*repmat(G_temp,1,L);
                    Y_temp = Y_temp+TG*V(:,:,i_Num);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update G
        Gn = permute(G(i,:,:),[2 3 1]);
        G(i,:,:) = max(Gn.*sqrt(permute(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N).*X_sp.*Y_sp.^(-2),2),[3 1 2])./permute(sum(repmat(T(:,:,i)*V(:,:,i),1,1,N)./Y_sp,2),[3 1 2])),0.01);
        for i_N = 1:N
            Y_temp = zeros(K_m,frame_N);
                for i_Num = 1:Num
                    G_temp = permute(G(i_Num,i_N,:),[3 1 2]);
                    TG = T(:,:,i_Num).*repmat(G_temp,1,L);
                    Y_temp = Y_temp+TG*V(:,:,i_Num);
                end
            Y_sp(:,:,i_N) = Y_temp;
        end
    end
    %%%%% IP of AuxIVA
    dlw = 0;
    for i = 1:K_m
        X_f = permute(X_Norm(i,:,:),[3 2 1]);
        dlw = dlw +log(abs(det(Q(:,:,i)))+epsi);
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
    Obj = ((sum(sum(sum(log(Y_sp+epsi))))+sum(sum(sum(X_sp./Y_sp))))/frame_N-2*dlw)/(Num*K_m);
    dObj = pObj-Obj;
    pObj = Obj;
    A(iteration,:) = [Obj,abs(dObj)/abs(Obj)];
%     if(abs(dObj)/abs(Obj)<theta)
%        break;
%     end
    %%%%%%%% Adjust the scales
    for i = 1:K_m
        for i_N = 1:N
            miu = Q(i_N,:,i)*Q(i_N,:,i)';
            Q(i_N,:,i) = Q(i_N,:,i)/sqrt(miu);
            G(:,i_N,i) = G(:,i_N,i)/miu;
        end
        for i_N = 1:Num
            phi = sum(G(i_N,:,i));
            G(i_N,:,i) = G(i_N,:,i)/phi;
            T(i,:,i_N) = T(i,:,i_N)*phi;
        end
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
            lamdag(i_N,:) = lamda(i,i_T,i_N)*G(i_N,:,i);
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