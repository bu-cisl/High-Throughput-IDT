%% Main_ID:
% This is the Main file to implement the transmission-mode intensity
% diffraction tomography(IDT) reconstruction algorithm.  The transfer
% functions (phase & amplitude) are derived based on the 1st Born
% approximation and the weak object approximation.
% Ref:
% L. Ruilong, et al. "High-throughput intensity diffraction tomography
% with a computational microscope." Biomed Opt Express (2018)

clc; clear;

%% measured data 
file_path = 'LED_algera_data';

%% add directory path of the functions
addpath(['.\Utilities']);

% if using gpu
gpu = 0;
% if performing lightfield refocus
lightfield = 0;
% if loading the data
loaddata = 1;

%% define commonly used functionsfr
Ft2 = @(x) fftshift(ifft2(ifftshift(x)));
% % % Define Fourier operators
F2 = @(x) fftshift(fft2(ifftshift(x)));
% % Define shift operator
% S = @(f,H) real(Ft(F(f).*H));
F = @(x) fft2(x);
Ft = @(x) ifft2(x);


%% experimental parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical aperture of the objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NA = 0.25;
lambda = 0.630; %in um
% maximum spatial frequency set by NA
um_m = NA/lambda;
% system resolution based on the NA
dx0 = 1/um_m/2;


%% parameters in real space
mag = 10.07; % microscope magnification
%dx0 = 6.5; % pixel size of camera
dpix_c = 6.5; %6.5um pixel size on the sensor plane
% effective image pixel size on the object plane
dpix_m = dpix_c/mag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% # of pixels at the output image patch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np = [1024,1024];
%Np = [1500,1000];
Nx = Np(1);
Ny = Np(2);

% FoV in the object space
FoV = Np*dpix_m;

dx = dpix_c/mag;
dy = dx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up image/object corrdinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncent = [1024,1024];
nstart = ncent-Np/2;
% start pixel of the image patch
img_center = (nstart-ncent+Np/2)*dpix_m;


%% index of surrounding media
n0 = 1.33;
%% parameters in Fourier space (DFT)
umax = 1/2/dpix_m;
dv = 1/dpix_m/Np(1); du = 1/dpix_m/Np(2);
u = [-umax:du:umax-du];
v = [-umax:dv:umax-dv];
[uu,vv] = meshgrid(u,v);
%% compute fz
eta = real(sqrt(1/lambda^2-(uu.^2+vv.^2)));

%% LED array geometries and derived quantities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spacing between neighboring LEDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ds_led = 4e3; %4mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from the LED to the object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_led = 79e3; % um

% size of each LED
s_led = 136; % um

% NA of single LED
NA_led = s_led/z_led;

% correponding spatial frequency is
fs_led = NA_led/lambda;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up LED coordinates
% h: horizontal, v: vertical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lit_cenv = 16;
lit_cenh = 16;
vled = [0:31]-lit_cenv;
hled = [0:31]-lit_cenh;

%% calculate dia_led
% total number of LEDs used in the experiment
[uu_led,vv_led] = meshgrid(hled,vled);
rrled = sqrt(uu_led.^2+vv_led.^2);

dd = sqrt((-uu_led*ds_led-img_center(1)).^2+(-vv_led*ds_led-img_center(2)).^2+z_led.^2);
sin_thetah = (-uu_led*ds_led-img_center(1))./dd; % axis cliped
sin_thetav = -(-vv_led*ds_led-img_center(2))./dd; % axis cliped

tan_thetah = (-uu_led*ds_led-img_center(1))./z_led;% axis cliped
tan_thetav = -(-vv_led*ds_led-img_center(2))./z_led;% axis cliped

%% determine brightfield (BF) or darkfield (DF)
% Model only works for BF
illumination_na = sqrt(sin_thetav.^2+sin_thetah.^2);
% corresponding spatial freq for each LEDs
for r0 = [5:0.01:25]
    s = rrled < r0;
    if (sum(s(:)) == sum(illumination_na(:)<=NA))
        dia_led = 2*r0;       
    end
end

LitCoord = rrled<dia_led/2;
% total number of LEDs used in the experiment

Nled = sum(LitCoord(:));
fvled = sin_thetav/lambda;
fuled = sin_thetah/lambda;
% spatial freq index for each plane wave relative to the center
idx_u = round(fuled/du);
idx_v = round(fvled/du);

illumination_na_used = illumination_na(LitCoord);
uu_prime_used = fuled(LitCoord);
vv_prime_used = fvled(LitCoord);

% number of brightfield image
NBF = sum(illumination_na_used<NA);
idx_BF = find(illumination_na_used<NA);
% darkfield images indices
NDF = sum(illumination_na_used>NA);
idx_DF = find(illumination_na_used>NA);


%% load only BF images
uu_BF = zeros(1,NBF);
vv_BF = zeros(1,NBF);

if loaddata
    I = zeros(Nx,Ny,length(idx_BF));
end
imgpath = sprintf('./%s/',file_path);
imglist = dir([imgpath,'*tif']);

%%%%%%%%%
y_idx = zeros(nnz(LitCoord),1);
x_idx = zeros(nnz(LitCoord),1);

count = 1;
for i = 1:32
    for j = 1:32
        if LitCoord(j,i) == 1
            img_array(:,:,count) = imread(sprintf('./%s\\img_y%d_x%d.tif',file_path,j,i));
            count = count+1;
        end
    end
end
%%%%%%%%%
n = 1;
for i = 1:32
    for j = 1:32
        if LitCoord(j,i) == 1
            nn = idx_BF(n);
            if loaddata
                
                fn = [imgpath,imglist(nn).name];
                tp = double(imread(sprintf('./%s\\img_y%d_x%d.tif',file_path,j,i)));
                I(:,:,n) = BkgndRemoval(tp,gpu);
            end
            uu_BF(n) = uu_prime_used(nn);
            vv_BF(n) = vv_prime_used(nn);
            n = n + 1;
        end
    end
end

clear tp dd idx_u idx_v illumination_na imglist rrled s img_array

%% Coordinates of the object / image

x = -(Np(2)-mod(Np(2),2))/2:1:(Np(2)-mod(Np(2),2))/2-(mod(Np(2),2)==0);
x = dpix_m*x;
y = -(Np(1)-mod(Np(1),2))/2:1:(Np(1)-mod(Np(1),2))/2-(mod(Np(1),2)==0);
y = dpix_m*y;
[xx,yy] = meshgrid(x,y);

%%
% reconstruction z-range
dz = 5;
z = [-20:dz:100]*n0;
Nz = length(z);
%% lightfield refocus: an intuitive way of visualizing the depth information
%% contained in the data
if lightfield
    showlf = 1;
    
    tan_thetah_used = tan_thetah(LitCoord);
    tan_thetav_used = tan_thetav(LitCoord);
    tan_theta = [tan_thetav_used(idx_BF),tan_thetah_used(idx_BF)];
    
    opts.method = 0; % use spatial domain method for fast shifting
    if opts.method == 0
        opts.Nobj = Np;
        opts.dx = dpix_m;
        opts.x = xx; opts.y = yy;
    else
        opts.F = F2; opts.Ft = Ft2;
        opts.u = uu; opts.v = vv;
    end
    I_lf = LF_refocus( I, tan_theta, z/n0, opts);
    if showlf
        for m = 1:Nz
            f1 = figure(1);
            imagesc(I_lf(:,:,m)); axis image; axis off;
            colormap gray;
            export_fig(f1,[imgpath(3:5),'_LF_z',num2str(round(z(m)/n0)),'.png'],'-m2');
        end
    end
    clear I_lf
end

%% the transfer function based IDT reconstruction starts from here
%% set up source function
eta_prime = real(sqrt(1/lambda^2-(uu_BF.^2+vv_BF.^2)));

% Generate source function Sf(u',v')
% size of each LED in the spatial frequency domain is
sfs_led = fs_led / du; % in diameter

Source = zeros(Nx,Ny,NBF);
Source_flip = zeros(Nx,Ny,NBF);

for m = 1:NBF
    % original source
    Source(:,:,m) = SourceFinitePtComp(uu_BF(m), vv_BF(m), uu, vv, sfs_led);
    % flipped source for transfer function computation
    Source_flip(:,:,m) = SourceFinitePtComp(-uu_BF(m), -vv_BF(m), uu, vv, sfs_led);
    %figure(81);imagesc(u,v,Sf(:,:,m));colormap jet;axis square;drawnow;
end


%% setup pupil function
% define pupil function P(u,v)
Pupil = sqrt(uu.^2+vv.^2)<=NA/lambda;
Pupil = double(Pupil);
% flipped pupil for transfer function computation
Pupil_flip = Pupil;



%%
% %% section 1: compute transfer functions at all z'
% % compute Hreal(u,v;z';u')

Hreal = zeros(Nx,Ny,Nz,NBF); % init
Himag = zeros(Nx,Ny,Nz,NBF); % init

w = waitbar(0,'Transfer Function Calculating...');
perc = 0;

savetransfunc = 0;
for m = 1:NBF
    S = Source(:,:,m);
    Sf = Source_flip(:,:,m);
    
    [Himag(:,:,:,m),Hreal(:,:,:,m)] = ...
        TransferFunction_3D_uvz3(lambda,S,Sf,Pupil,Pupil_flip,...
        z,dz,uu,vv,dx,n0);
    
    
    if savetransfunc
        %% plots to save
        % in-focus xy, Himag
        f1 = figure(1);imagesc(real(fftshift(Himag(:,:,round(Nz/2),m)))),colorbar;colormap jet;
        axis image; caxis([cbmin2,cbmax2]);axis off; drawnow;
        export_fig(f1,['samp_f1_',num2str(m),'.png'],'-m2');
        
        
        f4 = figure(4);imagesc(imag(fftshift(Hreal(:,:,round(Nz/2),m)))),colorbar;colormap jet;
        axis image; caxis([cbmin1,cbmax1]);axis off; drawnow;
        export_fig(f4,['samp_f4_',num2str(m),'.png'],'-m2');
        
    end
    %% display waitbar
    if mod(m,round((NBF-1)/10))==0
        perc = perc+10;
        waitbar(perc/100,w,sprintf('Transfer Function Calculating...%d%%',perc))
    end
    
end
close(w);


%% clean up memory:
clear Pupil Pupil_flip S Sf Source Source_flip sin_thetah sin_thetav tan_thetah tan_thetav
clear tan_thetah_used tan_thetav_used u v uu vv uu_BF vv_BF uu_led vv_led
clear uu_prime_used vv_prime_used xx yy x y fuled fvled idx_BF idx_DF illumination_na_used
clear eta eta_prime uled vled vled hled
%% Tikhonov regularization starts here

%% construct inversion transfer functions

HIreal = 0;
HIimag = 0;
sum_Himag = 0;
sum_Hreal = 0;
conj_part1 = 0;
conj_part2 = 0;

w = waitbar(0,'Calculating effective transfer function for Tikhonov Regularization...');
perc = 0;

for i = 1:NBF
    Itmp = I(:,:,i);
    Ihat_tmp = F(Itmp);
    
    HIreal  =  HIreal + conj(Hreal(:,:,:,i)) .* repmat(Ihat_tmp,1,1,Nz);
    HIimag  =  HIimag + conj(Himag(:,:,:,i)) .* repmat(Ihat_tmp,1,1,Nz);
    sum_Himag = sum_Himag + abs(Himag(:,:,:,i)).^2;
    sum_Hreal = sum_Hreal + abs(Hreal(:,:,:,i)).^2;
    conj_part1 = conj_part1 + Hreal(:,:,:,i).*conj(Himag(:,:,:,i));
    conj_part2 = conj_part2 + conj(Hreal(:,:,:,i)).*Himag(:,:,:,i);
    
    %% display waitbar
    if mod(i,round((NBF-1)/10))==0
        perc = perc+10;
        waitbar(perc/100,w,...
            sprintf('Calculating effective transfer function for Tikhonov Regularization...%d%%',perc))
    end
    
end
close(w);

%%
mu1 = 1 * 10.^(4);
mu2 = 1 * 10.^(4);

temp = (sum_Hreal + mu2).*(sum_Himag + mu1)- conj_part1.*conj_part2;

% dimension of V_im: fx, fy, z, mu
V_im = ((sum_Hreal + mu2).* HIimag - conj_part1.* HIreal) ./ temp;
V_re = ((sum_Himag + mu1).* HIreal - conj_part2.* HIimag) ./ temp;

v_im = real(Ft(V_im));
v_re = real(Ft(V_re));

n_re = sqrt(((n0.^2 + v_re) + sqrt((n0^2 + v_re).^2 + v_im.^2)) / 2);
n_im = v_im ./ n_re ./ 2;
cbmax1 = max([n_re(:)]);
cbmin1 = min([n_re(:)]);

cbmax2 = max([n_im(:)]);
cbmin2 = min([n_im(:)]);

phase_dir = sprintf('%s/phase_result',file_path);
absorp_dir = sprintf('%s/absorb_result',file_path);
mkdir(phase_dir);
mkdir(absorp_dir);
% save phase
for k = 1: Nz
    figure(1);imagesc(n_re(:,:,k));colormap gray;axis off;axis image;
    saveas(gcf,sprintf('%s/ph_slice%d.png',phase_dir,k));
end
% save absorption
for k = 1: Nz
    figure(2);imagesc(n_im(:,:,k));colormap gray;axis off;axis image;
    saveas(gcf,sprintf('%s/ab_slice%d.png',absorp_dir,k));
end



%% save results - n_re and n_im are the reconstructed object's ref-index
savefn = sprintf('algera_result_%s.mat',file_path);
save(savefn, 'NA', 'lambda', 'dpix_c', 'mag', 'n0', 'z_led', 's_led', ...
    'lit_cenh', 'lit_cenv', 'dia_led', 'n_re', 'n_im','v_im','v_re');