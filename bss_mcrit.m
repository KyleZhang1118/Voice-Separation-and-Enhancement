function [SDR,ISR,SIR,SAR]=bss_mcrit(s_true,e_spat,e_interf,e_artif)

% BSS_MCRIT Computation of some separation performance criteria given the
% decomposition of an estimated source image into four components
% representing respectively the true source image, spatial (or filtering)
% distortion, interference and artifacts.
%
% [SDR,ISR,SIR,SAR]=bss_mcrit(s_true,e_spat,e_interf,e_artif)
%
% Inputs:
% s_true: I x T matrix containing the true source image (one row per channel)
% e_spat: I x T matrix containing the spatial (or filtering) distortion component
% e_interf: I x T matrix containing the interference component
% e_artif: I x T matrix containing the artifacts component
%
% Outputs:
% SDR: Signal to Distortion Ratio
% ISR: source Image to Spatial distortion Ratio
% SIR: Source to Interference Ratio
% SAR: Sources to Artifacts Ratio

%%% Errors %%%
if nargin<4, error('Not enough input arguments.'); end
[It,Tt]=size(s_true);
[Is,Ts]=size(e_spat);
[Ii,Ti]=size(e_interf);
[Ia,Ta]=size(e_artif);
if ~((It==Is)&&(It==Ii)&&(It==Ia)), error('All the components must have the same number of channels.'); end
if ~((Tt==Ts)&&(Tt==Ti)&&(Tt==Ta)), error('All the components must have the same duration.'); end

%%% Energy ratios %%%
% SDR
SDR=10*log10(sum(sum(s_true.^2))/sum(sum((e_spat+e_interf+e_artif).^2)));
% ISR
ISR=10*log10(sum(sum(s_true.^2))/sum(sum(e_spat.^2)));
% SIR
SIR=10*log10(sum(sum((s_true+e_spat).^2))/sum(sum(e_interf.^2)));
% SAR
SAR=10*log10(sum(sum((s_true+e_spat+e_interf).^2))/sum(sum(e_artif.^2)));

return;