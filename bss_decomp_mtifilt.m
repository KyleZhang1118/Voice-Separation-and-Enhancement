function [s_true,e_spat,e_interf,e_artif]=bss_decomp_mtifilt(se,S,j,L)

% BSS_DECOMP_MTIFILT Decomposition of an estimated source image into four
% components representing respectively the true source image, spatial (or
% filtering) distortion, interference and artifacts, derived from the true
% source images using multichannel time-invariant filters.
%
% [s_true,e_spat,e_interf,e_artif]=bss_decomp_mtifilt(se,S,j,L)
%
% Inputs:
% se: I x T matrix containing the estimated source image (one row per channel)
% S: I x T x J matrix containing the true source images
% j: source index corresponding to the estimated source image in S
% L: length of the multichannel time-invariant filters in samples
%
% Outputs:
% s_true: I x T matrix containing the true source image, i.e. S(:,:,j) (one row per channel)
% e_spat: I x T matrix containing the spatial (or filtering) distortion component
% e_interf: I x T matrix containing the interference component
% e_artif: I x T matrix containing the artifacts component

%%% Errors %%%
if nargin<4, error('Not enough input arguments.'); end
[Ie,Te]=size(se);
[I,T,J]=size(S);
if I~=Ie, error('The number of channels of the true source images and the estimated source image must be equal.'); end
if T~=Te, error('The duration of the true source images and the estimated source image must be equal.'); end

%%% Decomposition %%%
% True source image
s_true=[S(:,:,j),zeros(I,L-1)];
% Spatial (or filtering) distortion
e_spat=project(se,S(:,:,j),L)-s_true;
% Interference
e_interf=project(se,S,L)-s_true-e_spat;
% Artifacts
e_artif=[se,zeros(I,L-1)]-s_true-e_spat-e_interf;

return;



function sproj=project(se,S,L)

% Least-squares projection of each channel of se on the subspace spanned by
% delayed versions of the channels of S, with delays between 0 and L-1

[I,T,J]=size(S);
S=reshape(permute(S,[1 3 2]),I*J,T);

%%% Computing coefficients of least squares problem via FFT %%%
% Zero padding and FFT of input data
S=[S,zeros(I*J,L-1)];
se=[se,zeros(I,L-1)];
N=2^nextpow2(T+L-1);
Sf=fft(S,N,2);
sef=fft(se,N,2);
% Inner products between delayed versions of S
G=zeros(I*J*L);
for k1=0:I*J-1
    for k2=0:k1
        SSf=Sf(k1+1,:).*conj(Sf(k2+1,:));
        SSf=real(ifft(SSf));
        SS=toeplitz(SSf([1 N:-1:N-L+2]),SSf(1:L));
        G(k1*L+1:k1*L+L,k2*L+1:k2*L+L)=SS;
        G(k2*L+1:k2*L+L,k1*L+1:k1*L+L)=SS.';
    end
end
% Inner products between se and delayed versions of S
D=zeros(I*J*L,I);
for k=0:I*J-1
    for i=1:I
        Ssef=Sf(k+1,:).*conj(sef(i,:));
        Ssef=real(ifft(Ssef,[],2));
        D(k*L+1:k*L+L,i)=Ssef(:,[1 N:-1:N-L+2]).';
    end
end

%%% Computing projection %%%
% Distortion filters
C=G\D;
C=reshape(C,L,I*J,I);
% Filtering
sproj=zeros(I,T+L-1);
for k=1:I*J
    for i=1:I
        sproj(i,:)=sproj(i,:)+fftfilt(C(:,k,i).',S(k,:));
    end
end

return;