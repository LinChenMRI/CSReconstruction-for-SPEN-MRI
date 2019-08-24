%% Random sampling and CS-based reconstruction algorithm for SPEN MRI
% Remove the aliasing artifacts via random sampling and CS-based reconstruction
% Lin Chen
% Xiamen University
% Email: chenlin@stu.xmu.edu.cn
% Aug. 5, 2015
% If you use this code, please cite the following papers:
% [1] L. Chen, L.J. Bao, J. Li, S.H. Cai, C.B. Cai, Z. Chen, An aliasing artifacts reducing approach with random undersampling for spatiotemporally encoded single-shot MRI, J. Magn. Reson., 237 (2013) 115-124.

clear all; close all; clc
addpath(strcat(pwd,'/Toolbox'));
%%
% * Create the time stamp*
date_now=datestr(now);
date_hour=num2str(hour(date_now));
date_minute=num2str(minute(date_now));
date_second=num2str(second(date_now));
time=['_',date_hour,'_',date_minute,'_',date_second];

%%
% * Choose the sampling pattern. 
% 0 for Equispaced sampling and 1 for random sampling
param.samplingpattern = 1;
param.fft_number = 256; % The digital resolution of SR image
param.SVDpreprocessing = 0; % SVD is used to reduce the condition number of linear equations to make the algorithm more robust
% * Input the directory of the FID file* 
if param.samplingpattern ==0
    param.fid_dir = 'Data\Equispaced sampling'; % The directory of FID file
else 
    param.fid_dir = 'Data\Random sampling';
end
Procpar = readprocpar(param.fid_dir);   % load the experiment parameters
%% Get the experiment parameters

% * Some common paramters* 
param.gama = 42.58e-6;
param.segment = 1;  % The number of segments used in the experiment

% * Extract the paramters which are needed in the reconstruction process* 
param.Lpe = Procpar.lpe.Values*1e-2; % The FOV of phase encoded dimension
param.pw = Procpar.p1_bw.Values; % The bandwidth of the chirp pulse
param.Texc = Procpar.pw1.Values; % The duration of the chirp pulse
param.nphase = Procpar.nphase.Values; % The number of points in the phase encoded dimension
param.nread = Procpar.nread.Values/2; % The number of points in the frequency encoded dimension
param.chirp = 1; % 1 for 90 chirp and 2 for 180 chirp. 
param.sign = 1; % choose 1 or -1 according to the sweep direction of chirp

% * Calculate the other paramters* 
param.Gexc = param.sign*param.pw/(param.gama*param.Lpe);
param.gama = param.gama*2*pi; % translate the gama into rid form
%%
% *Load the original data in the SPEN sampling domain*
[RE,IM,NP,NB,NT,HDR]=load_fid(param.fid_dir,param.segment); % load fid data,choose the wanting segment
fid=RE+1i*IM;   % translate into complex form
fid=fid(1+param.nread*2:end,:); % discard the front useless data
zerofillnum=(param.fft_number-param.nread)/2;% the number of zero-filling points 

for n = 1:1:NT
    fid_temp = reshape(fid(:,n),param.nread,param.nphase).'; % Translate into 2D
    fid_temp(2:2:end,:) =  fliplr(fid_temp(2:2:end,:)); % Rearrange the even lines
    fid_temp_zerofill = [zeros(param.nphase,zerofillnum),fid_temp,zeros(param.nphase,zerofillnum)]; % Zero-fill the fid
    fid_temp_fftshift  = fftshift(fid_temp_zerofill,2);
    st_temp_fftshift = fft(fid_temp_fftshift,[],2);
    st_temp = fftshift(st_temp_fftshift,2);
    param.st(:,:,n) = st_temp;
end
figure(1);imshow(abs(param.st),[]);title('Blurred image');drawnow
%%
% *Create the directory for saving the results*
param.savedatadir = [param.fid_dir,'\result'];
if exist(param.savedatadir,'file')==0
    mkdir(param.savedatadir);
end

%%
% * SR reconstruction for the SPEN data* 
if param.samplingpattern ==0
    param.SPP_position = linspace(param.Lpe/2,-param.Lpe/2,param.nphase).'; % the position of stationary phase point
else
    load([param.fid_dir,'\samplingtrajectory']);
    param.samplingtrajectory = samplingtrajectory;
    param.SPP_position = -param.Lpe/2 + cumsum(param.samplingtrajectory)*param.Lpe/sum(param.samplingtrajectory);
end
param.interpolation_poistion = linspace(param.Lpe/2,-param.Lpe/2,param.fft_number).'; % the position of interpolated points

param.additional_term = param.gama*param.chirp*param.Gexc*param.Texc/param.Lpe;
param.addtional_phase = exp(1i*param.additional_term*param.SPP_position.^2/2);

for n = 1:1:size(param.st,3)
    param.st(:,:,n) = param.st(:,:,n).*repmat(param.addtional_phase,[1,param.fft_number]);
    for m = 1:1:param.fft_number
        interpolation_temp(:,m,n) = interp1(param.SPP_position,param.st(:,m,n),param.interpolation_poistion);
    end
end
param.st = interpolation_temp;
param.st(isnan(abs(param.st(:,:,:))))=0;

param.P_matrix = zeros(param.fft_number,param.fft_number);
for n = 1:1:param.fft_number
    for m = 1:1:param.fft_number
        param.P_matrix(n,m)=exp(1i*param.additional_term*(param.interpolation_poistion(m)-param.interpolation_poistion(n))^2/2);
    end
end

if param.SVDpreprocessing ==1   % whether utilize the SVD preprocessing
    [U,S,V] = svd(param.P_matrix);
    param.P_matrix = U*eye(size(S))*V';
end
param.P_matrix = Pmatrix(param.P_matrix);
result_without_CS = param.P_matrix'*param.st;
result_without_CS = imrotate(result_without_CS,90);
figure(2);imshow(abs(result_without_CS),[]);title('SR result without CS');drawnow
imwrite(abs(result_without_CS)/max(abs(result_without_CS(:))),[param.savedatadir,'\result_without_CS.tiff'],'tiff')

CSparam=init;
CSparam.FT=Pmatrix(param.P_matrix);
CSparam.data=param.st/norm(param.st);
CSparam.Itnlim=20;
CSparam.TV=TVOP;
CSparam.XFM=1;
% CSparam.XFM=Wavelet('Daubechies',6,2);
CSparam.xfmWeight=0.01; % The regularization parameters may need to be fine-tuned
CSparam.TVWeight=0.01;
x=zeros([param.fft_number,param.fft_number]);

for n=1:1:5
    x=fnlCg(x,CSparam);
    figure(100), imshow(imrotate(abs(x),90),[]), drawnow
end
result_with_CS = x;
result_with_CS = imrotate(result_with_CS,90);
figure(3);imshow(abs(result_with_CS),[]);title('SR result with CS');
imwrite(abs(result_with_CS)/max(abs(result_with_CS(:))),[param.savedatadir,'\result_with_CS.tiff'],'tiff')