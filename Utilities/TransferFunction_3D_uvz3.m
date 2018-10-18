% THIS IS ACTUALLY uvz4
function [ImagTransFunc,RealTransFunc] = ...
    TransferFunction_3D_uvz3(lambda,S,Sf,P,Pf,z,dz,Fx,Fy,ps,RI)
% illumination-dependent normalized Weak Object Transfer Function
%
% Calculate the illumination-dependent transfer functions at each depth slices
% for both real and imaginary parts of the permittivity contrast
%
%
% Inputs:
%       lambda = wavelength of illumination
%       S = source function
%       Sf = flipped source function
%       P  = pupil function
%       Pf = flipped pupil function
%       z = position of each layer
%       dz = sample size along the axial direction
%       Fx = x coordinate in Fourier space
%       Fy = y coordinate in Fourier space
%       ps = pixel size in real space
%       RI = refactiv index of the surrounding medium
%       padsize = zeropadding size in the real space
%
% Outputs:
%       ImaginaryTransferFunction = weak absorption transfer function
%       RealTransferFunction     = weak phase transfer function
%
% adapted from Michael Chen, Sep 2015
% define function handle of Fouirer transform and inverse
%
% Last modified by Lei Tian (lei_tian@alum.mit.edu), 07/10/2017
% fixed an important const missing in the code


F = @(x) fft2(x);
IF = @(x) ifft2(x);

% number of pixels
[M,N] = size(Fx);

% constants in the transfer functions
%% Lei's note: Lei found a few factors are missing in the code
k0 = 2*pi/lambda;
const = (1/2)*k0^2*dz;  

%% modified by Lei Tian on 7/6/17
% take real part to remove the evanescent components
G = real(1./(2*pi*sqrt((RI/lambda)^2-(Fx.^2+Fy.^2))));%%


% init
ImagTransFunc = zeros(M,N,length(z));
RealTransFunc = zeros(M,N,length(z));

for j = 1:length(z)
    
    %% Lei Tian modified on 7/10/2017
    prop_phase = exp(1i*2*pi*z(j).*real(sqrt((RI/lambda)^2-(Fx.^2+Fy.^2))));

    FPSfph_cFPphG = F(S.*conj(Pf).*prop_phase).*F(P.*conj(prop_phase).*G);
    FPSfph_cFPphG = ifftshift(fftshift(FPSfph_cFPphG));
    ImagTransFunc(:,:,j) = IF(real(FPSfph_cFPphG));
    RealTransFunc(:,:,j) = IF(1i*imag(FPSfph_cFPphG));
    
    
end

%% compute the DC (I1) term for normalization
DC = sum(sum(abs(P).^2.*Sf));

%% compute normalized transfer functions
ImagTransFunc = -const*ImagTransFunc./DC;
RealTransFunc = 1i*const*RealTransFunc./DC;

% remove points due to numerical errors
RealTransFunc(abs(RealTransFunc)<1e-3) = 0;
ImagTransFunc(abs(ImagTransFunc)<1e-3) = 0;


end

