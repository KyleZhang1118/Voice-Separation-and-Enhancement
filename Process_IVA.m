function [Y,W,SetupStruc] = Process_IVA(s,Transfer,SetupStruc)
K = SetupStruc.IVA.K;
hop = SetupStruc.IVA.hop;
win = hanning(K,'periodic');
win = win/sqrt(sum(win(1:hop:K).^2));
SetupStruc.IVA.win = win;  % Preserve 'win' in 'SetupStruc'
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
theta = 10^-4;
for i = 1:K_m
    X_f = permute(X(i,:,:),[3 2 1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize W by PCA
    if(i>1)
        Steer = permute(Transfer(i,:,:),[2 3 1]);
    else 
        Steer = ones(N,Num);
    end
    [E,D] = PCA(X_f,1,Num);
    V = sqrt(D)\E';
    V_sp(:,:,i) = V;
    X_sp(:,:,i) = V*X_f;
    Ori = V*Steer;
    if rcond(Ori)<theta
        Ori = Ori+eye(Num)*min(diag(Ori));
    end
    %%%%%%%%%%%% Adjust amplitude of 'w'
    W_o = inv(Ori);
    y_f = W_o*V*X_f;
    norm = max(abs(y_f),[],2);
    if(norm>10)
        norm = repmat(norm,1,Num);
        W_o = W_o./norm;
        y_f = W_o*V*X_f;
    end
    W_IVA(:,:,i) = W_o;
    Y_f(:,:,i) = y_f;    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% IVA iterations
step_size = 0.1;
max_iteration = 1000;
Y_k = zeros(Num,frame_N);
K_sp = sqrt(K_m);
i_sign = zeros(1,K_m)+1;
state_sign = zeros(1,K_m);
A = zeros(1001,K/2)-1; %%%% Show the decrease of non-linear correlation, IVA max iterations 1000
for iteration = 1:max_iteration
    for i = 1:Num
        y_temp = permute(Y_f(i,:,:),[3 2 1]);
        Y_k(i,:) = sqrt(sum(abs(y_temp)));
    end
    for i = 1:K_m
        W = W_IVA(:,:,i);
        X_f = X_sp(:,:,i);
        y_f = Y_f(:,:,i);
        y_fun = K_sp*y_f./Y_k;
        core = eye(size(W))-y_fun*y_f'/frame_N;
        if(det(core)<1e-4)
            state_sign(i) = 1;
            continue;
        end
        state_sign(i) = 0;
        W = W+step_size*core*W;
        W_IVA(:,:,i) = W;
        Y_f(:,:,i) = W*X_f;
        A(i_sign(i),i) = abs(det(core));
        i_sign(i) = i_sign(i)+1;
    end
    if(~ismember(0,state_sign))
        break;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Post processing
W = zeros(Num,N,K_m);
Y_f(:,:,1) = zeros(Num,frame_N);
for i = 2:K_m
    Steer = permute(Transfer(i,:,:),[2 3 1]);
    I_pro = W_IVA(:,:,i)*V_sp(:,:,i)*Steer;
    for ii = 1:Num
        Y_f(ii,:,i) = Y_f(ii,:,i)/I_pro(ii,ii);
        W_IVA(ii,:,i) = W_IVA(ii,:,i)/I_pro(ii,ii);
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