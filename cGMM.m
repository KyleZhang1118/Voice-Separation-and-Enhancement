function [lamda,R,Q] = cGMM(X,Transfer,iters,R_initial)
if(~exist('iters','var'))
    iters = 20;
end
K_m = size(X,1);
frame_N = size(X,2);
N = size(X,3);
Num = size(Transfer,3);
Num_v = Num+1;
theta = 1e-4;
lamda = zeros(frame_N,Num_v,K_m);
phi = zeros(frame_N,Num_v);
R_X = zeros(N,N,frame_N);
Q = zeros(iters,K_m);
p = zeros(frame_N,Num_v,K_m);
const = (2*pi)^(N/2);
R = zeros(N,N,Num_v,K_m);
for i = 1:K_m
    if(~exist('R_intial','var'))
        RTF = permute(Transfer(i,:,:),[2 3 1]);
        for j = 1:Num
            R(:,:,j,i) = RTF(:,j)*RTF(:,j)';
        end
    else
        for j = 1:Num
            R(:,:,j,i) = R_initial(:,:,i,j);
        end
    end
    if(Num_v == Num+1)
        R(:,:,Num_v,i) = eye(N);
    end
end
PrintLoopPCw('   Starting iterations of cGMM. ');
for i = 1:K_m
%     if(i==13)
%         tt =1;
%     end
    X_i = permute(X(i,:,:),[3 2 1]);
    R_i = R(:,:,:,i);
    lamda_i = zeros(frame_N,Num_v);
    p_i = zeros(frame_N,Num_v);
    alpha = ones(1,Num_v)/(Num_v);
    for t = 1:frame_N
        R_X(:,:,t) = X_i(:,t)*X_i(:,t)';
    end
    for iters_N = 1:iters
        for j = 1:Num_v
            if(rcond(R_i(:,:,j))<theta)
                R_i(:,:,j) = R_i(:,:,j)+theta*eye(N)*max(diag(R_i(:,:,j)));
            end
            for t = 1:frame_N
%                 phi(t,j) = sum(diag(R_X(:,:,t)))/sum(diag(R_i(:,:,j)));
                phi(t,j) = sum(diag(R_X(:,:,t)/R_i(:,:,j)))/N;
                temp = R_i(:,:,j)*phi(t,j);
%                 p_i(t,j) = abs(exp(-X_i(:,t)'/temp*X_i(:,t)/2)/(sqrt(det(temp))*const));
%                 const is a constant which doesnt affect the ratio between the components
                p_i(t,j) = abs(exp(-X_i(:,t)'/temp*X_i(:,t)/2)/sqrt(det(temp)));
            end
        end
        for t = 1:frame_N
            lamda_i(t,:) = alpha.*p_i(t,:);
            lamda_i(t,:) = lamda_i(t,:)/sum(lamda_i(t,:));
        end
        Q(iters_N,i) = sum(sum(lamda_i.*log(repmat(alpha,frame_N,1).*p_i)));
        for j = 1:Num_v
            R_temp = zeros(N);
            for t = 1:frame_N
                R_temp = R_temp+lamda_i(t,j)*R_X(:,:,t)/phi(t,j);
            end
            R_i(:,:,j) = R_temp/(sum(lamda_i(:,j))+eps);
            alpha(j) = sum(lamda_i(:,j))/frame_N;
        end
    end
    p(:,:,i) = p_i;
    lamda(:,:,i) = lamda_i;
    R(:,:,:,i) = R_i;
    PrintLoopPCw(i,K_m);
end
for j = 1:Num_v
    figure
    lamda_j = permute(lamda(:,j,:),[3 1 2]);
    imagesc(flipud(lamda_j))
end
return














