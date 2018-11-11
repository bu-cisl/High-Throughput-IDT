% illumination-dependent Weak Object Transfer Function
% 
% Calculate the illumination-dependent transfer functions at each depth slices
% for both real and imaginary parts of the permittivity contrast
%
%
% Inputs:
%       lambda = wavelength of illumination
%       S = source shape
%       P  = pupil function
%       z = position of each layer
%       dz = sample size along the axial direction
%       Fx = x coordinate in Fourier space
%       Fy = y coordinate in Fourier space
%       ps = pixel size in real space
%       RI = refactiv index of the surrounding medium
%
% Outputs:
%       ImaginaryTransferFunction = weak absorption transfer function      
%       RealTransferFunction     = weak phase transfer function
%

function [ImagTransFunc,RealTransFunc] = ...
    TransferFunction_3D_uvz2(lambda,S,P,z,dz,Fx,Fy,ps,RI)
    
    % define function handle of Fouirer transform and inverse
    F = @(x) fft2(x);
    IF = @(x) ifft2(x);

    % sampling size the Fourier space
    dfx = Fx(1,2)-Fx(1,1);
    dfy = Fy(2,1)-Fy(1,1);
    
    % number of pixels
    [M,N] = size(Fx);
    
    % constants in the transfer functions
    %const = -1/lambda^2*dz;
    const = -1/lambda^2*dz;
    
%% include the source flip here. need to fix
%     Sf = padarray(S,[1,1],'post');
%     Sf = flip(Sf,1); Sf = flip(Sf,2);
%     Sf = Sf(1:end-1,1:end-1);
%     Sf = ifftshift(Sf);
    Sfftshift = ifftshift(S);
    %% modified by Lei Tian on 7/6/17
    % take real part to remove the evanescent components
    G = real(1./sqrt((RI/lambda)^2-(Fx.^2+Fy.^2)));%%
    
    w = waitbar(0,'Transfer Function Calculating...');
    perc = 0;
    
    % init
    ImagTransFunc = zeros(M,N,length(z));
    RealTransFunc = zeros(M,N,length(z));
    
    for j = 1:length(z)
        
        %% Lei Tian modified on 4/6/2017
        prop_phase = exp(1i*2*pi*z(j)*sqrt((RI/lambda)^2-(Fx.^2+Fy.^2)));
        prop_phase(abs(P)==0) = 0;  
        %FPSfph_cFPphG = F(Sf.*P.*prop_phase).*conj(F(P.*prop_phase.*G))*dfx*dfy;
        %FPSfph_cFPphG = F(Sf.*P.*prop_phase).*conj(F(P.*prop_phase.*G))*dfx*dfy;
        %% check the conj of P to see if they match Micheal's derivation
        FPSfph_cFPphG = F(Sfftshift.*conj(P).*prop_phase).*F(P.*conj(prop_phase).*G)*dfx*dfy;
        
        %% check the M*N term
        ImagTransFunc(:,:,j) = IF(real(FPSfph_cFPphG))*M*N*ps^2;
        RealTransFunc(:,:,j) = IF(1i*imag(FPSfph_cFPphG))*M*N*ps^2;        
    
        if mod(j,round((length(z)-1)/5))==0
           perc = perc+20;
           waitbar(perc/100,w,sprintf('Transfer Function Calculating...%d%%',perc))
        end

    end

%     ImaginaryTransferFunction = fft(ImaginaryTransferFunction,[],3)*(z(2)-z(1));
%     RealTransferFunction = fft(RealTransferFunction,[],3)*(z(2)-z(1));
%% DC is |P(kxi,kyi)|^2, need to correct this line
    %DC = sum(abs(P(:)).^2);
    DC = sum(sum(abs(P).^2.*Sfftshift));
    %DC = 1;
    %DC = abs(P(:,:)).^2;
    
    %% check if const is correct.
    ImagTransFunc = const*ImagTransFunc./DC;
    RealTransFunc = -1i*const*RealTransFunc./DC;

    close(w);

end

