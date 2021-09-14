# Voice-Separation-and-Enhancement
## Discription  
This program consists of several popular methods and its variants for speech separation and enhancement. The purposes of this program are to realize, test and compare methods quickly. The default model of microphone array is 6+1(peripheral+central) circular array. Test data are generated by ISM method[1,2] based on TIMIT database. Voicebox toolbox is required.  
All codes are written and updated in Matlab by Ke Zhang. If you find any bug or error, please contact me.(zhangke9311@163.com)  
The list of main methods:  
Beamforming:  
* DSB  
* MVDR  
* LCMV  
* maxSNR/GEVD  

Blind source separation(BSS):  
* ICA  
* IVA  
* AuxIVA  
* OverIVA  
* ILRMA  
* FastMNMF  

In general, methods in beamforming use the steering vector or other spatial information of sources to enhance the target speech, and BSS methods only use the number of sources except some cases for solving the permutation ambiguity.  
## Guides for users  
1. The main function is command.m in which you can set the number and angles of sound sources(0-45-315 degrees), and select the algorithms in the list you want to test(set the value behind the corresponding method to 1 for running, 0 for not). Simulation environment can be set in ISM_setup.m, such as T60 for reverberation(0, 0.3s, 0.6s, 0.9s support), configurations of microphone array and NoiseFlag for noise adding,etc.  
2. In Process.m, the steering vector of sound sources are calculated by the function 'Cal_transfer' and the length of window for fft is set. The plotting and generation of separated signals are controlled by 'sign_plot' and 'sign_write' in Process.m. The performance of methods are shown in the top right corner of the plotting figure where the column with 2 numbers consists of SIR improved and SIR output and the column with 3 numbers consists of SDR, SIR and SAR, respectively.   
3. If you want to test the method written by yourself, follow the steps:  
* step1: Add your method to the list in command.m  
* step2: Add the session of your method in Process.m  
* step3: Write the function of your method which often named as 'Process_name', the structure of which can refer to Process_DSB(a simplest method) which consists of the fundermental flow, such as the stft and the reconstruction of the separated speech from time-frequency domain to waveform in time domain. The output of the 'Process_name' should consist of 'Y' which is the separated speech, 'W' which is the filter or demixing matrix, and the structure variable 'SetupStruc' which consists of some essential parameters(noting the window of enframe should be reserved in structure variable named as 'SetupStruc.name' in 'Process_name').  
4. The room size used in ISM mehod is 6x4x3m. The microphone array is placed in the center of room and the distance between adjacent microphone in 6+1 microphone array is 4.35cm. The distance between the sound source and the center of microphone array is 1.5m. If you want to test other environments, you sould read and rewrite the readData.m to use the data generated by yourself. The calculation of steering vector in Cal_transefer.m should also be rewritten based on the configuration of the microphone array.  
5. CommandOnline.m and OnProcess.m are designed for online methods which process input signals by blocks, not whole. They are seldom used and not updated frequently.  

## Details of methods in realization 
### Beamforming  
#### DSB(Delay and sum)
* DSB is the simplest method for speech enhancement by compensating the phase difference between the first channel and others in time-frequency domain.
* DSB_Mask uses the steering vector of all sound sources, and is time-varying filter which use a post binary mask to only preserve the beam with the highest energy. The method is simple and effective to wipe off interference but lose the components of target speech. [3]  
#### MVDR and LCMV
* MVDR only uses the steering vector of target source, while LCMV uses all steering vectors to suppress the interference.  
* MVDR_ESB and LCMV_ESB use PCA to obtain the eigen subspace(whose dimension is the number of sound sources) of the covariance matrix, then do normal beamforming on the subspace matrix.  
* MVDR_AESB and LCMV_AESB are similar to MVDR_ESB, but the dimensions of subspace are determined by judgement of eigen values, not directly equalling to the number of sound sources.   
* MVDR_PCA is similar to MVDR_ESB, but without the diagonal loading on the subspace matrix.(This program is written to evaluate the sparsity of the speech components in time-frequency domain.)  
* MVDR_Search and LCMV_Search use the diagonal loading to produce the irreversible convariance matrix. The value of diagonal loading is searched by meeting the requirement of the WNG(White noise gain). [4]  
* MVDR_AESB_Search is a combination of MVDR_AESB and MVDR_Search, in which the value of diagonal loading is searched and put on the subspace matrix.  
#### maxSNR/GEVD  
* maxSNR uses the convariance matrices of all sound sources to calculate the maxSNR filters. These filters maybe can be considered as the best demixing matrix under the current configuration of microphone array in terms of SIR, but result in some loss of components of target speech. [5]  
* cGMM_maxSNR is similar to the maxSNR, while the convariance matrices are estimated by cGMM method. The intial of cGMM is calculated by the steering vector of sound sources to solve the ambiguity of permutation. [5,6]  
### Blind source separation(BSS)
#### ICA
* ICA_funda is the regular frequency-domain ICA realized on the information maximization approach combined with natural gradient. The permutation is regularized by calculating the correlation between the separated signal and the original signal of targer speech in frequency domain.  
* ICA_initial is similar to ICA_funda, while the dimixing matrix is initialized by the inverse matrix of the assumed mixing matrix which loads steering vectors. The post regularization for the ambiguity of permutation is not processed.  
* ICA_Sawada is similar to ICA_funda, while the post proceeding for the ambiguity of permutation is based on the Sawada method. [7]  
#### FastICA
* FastICA_HO_Sawada is frequency-domain FastICA based on the higher-order statistics, in which PCA has been used to reduce dimension and the ambiguity of permutation is sovled by the Sawada method. The separation vectors are iterated alternately and parallel orthogonalized.[8]  
#### IVA
* IVA is based on the minimization of the Kullback-Leibler divergence and the model of signal is spherically symmetric multivariate Laplace distribution. [9]  
* IVA_woDR is similar to IVA, while do not uses PCA to reduce the dimension of signals. The separated signals are selected by the energy.  
* AuxIVA uses the iterative projection(IP) algorithm to obtain the demixing matrix. [10]  
* OverIVA is designed for the over-determined situation, which does not use the PCA to reduce dimension. The dimixing matrix in OverIVA consists of two parts, one is searched by iterative projection(IP) to separate the target speech, another is an orthogonal subspace to help IP. OverIVA obtains a high SIR compared to AuxIVA which uses PCA, but losses the components of targer speech resulting in a low SAR. [11]  
#### ILRMA
* ILRMA is a combination of NMF and AuxIVA which uses the time-varying gaussian distribution and NMF to model the spectrogram of target sources. The improvement of ILRMA compared to AuxIVA is around 2dB in my test. [12]  
* ILRMA_woDR is similar to ILRMA, but does not use the PCA to reduce dimension. The purpose of ILRMA_woDR is to test the numerical stability of ILRMA. The separated signals are selected by the energy.  
* ILRMA_PF is similar to ILRMA, while the basis of NMF is shared across the signals and a partition function is used to compose the spectogram of target signals. [12]  
* OverILRMA is a over-determined version of ILRMA. Compared to OverIVA, NMF is used to model the spectrogram of target signals in the way which is similar in the ILRMA.  
#### FastMNMF
* FastMNMF1 and FastMNMF2 are the fast version of MNMF based on the jointly-diagonalizable(JD) full-rank spatial models. The spatial convariance matrices of target signals are composed by M(which commonly is the number of channels) rank-1 matrices. These methods are expected to work in the reverberation environments but not good in my test. The iteration seems not stable and is numerical sensitive. [13]  

## References
[1] Lehmann, E.A. Diffuse reverberation model for efficient image-source simulation of room impulse responses. 2009.  
[2] Available online: http://www.eric-lehmann.com/  
[3] Anatasios A. Capturing and reproducing spatial audio based on a circular microphone array. 2013.  
[4] Henry, C. Robust Adaptive Beamforming. 1987.  
[5] Ernst, W. Blind acoustic beamforming based on generalized eigenvalue decomposition. 2007.  
[6] Takuya, H. Online MVDR beamformer based on complex Gaussian mixture model with spatial prior for noise robust ASR. 2017.  
[7] Hiroshi, A. Robust and precise method for solving the permutation problem of frequency-domain blind source separation. 2004.  
[8] E. Bingham. A fast fixed-point algorithm for indenpendent component analysis of complex valued signals. 2000.  
[9] Taesu, K. Blind source separation exploiting higher-order frequency denpendencies. 2007.  
[10] Nobutaka, O. Stable and fast update rules for independent vector analysis based on auxiliary function technique. 2011.  
[11] Robin, S. Independent vector analysis with more microphones than sources. 2019.  
[12] Daichi, K. Determined blind source separation with independent low-rank matrix analysis. 2018.  
[13] Kouher, S. Fast multichannel nonnegative matrix factorization with directivity-aware jointly-diagonalizable spatial covariance matrices for blind source separation. 2020.  

Last edited in 9/14/2021
