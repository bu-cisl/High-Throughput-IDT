function [ Ibf, Imin, Imax ] = LF_refocus( I, tan_theta, z, opts)
%LF_REFOCUS performs lightfield refocus using different methods
%   method
%   Detailed explanation goes here
%
% refs:
%   [1] Lei Tian, Jingyan Wang, and Laura Waller,
%    "3D differential phase-contrast microscopy with computational
%    illumination using an LED array," Opt. Lett. 39, 1326-1329 (2014).
%
% last modified by Lei Tian (lei_tian@alum.mit.edu) 7/10/2017 @ Botston
% Univ.


Nslice = length(z);
[N1,N2,Nangle] = size(I);


if opts.method == 0 % use spatial domain
    lf_stack = I;
else
    % % Define Fourier operators
    F = opts.F;
    Ft = opts.Ft;
    Ift = F(I);
    Ift = padarray(Ift,[(opts.Nobj(1)-N1)/2,(opts.Nobj(2)-N2)/2]);
    lf_stack = Ift;
end

clear I Ift

%% define shift operators
if opts.method == 0 % use spatial domain
    lf_shift = @(lf, xs, ys) circshift(lf,...
        [-round(ys/opts.dx),-round(xs/opts.dx)]);
    
else % Fourier domain
    lf_shift = @(LF, xs, ys) ...
        real(Ft(LF.*exp(1i*2*pi*(opts.u*xs+opts.v*ys))));
end


Ibf = zeros(opts.Nobj(1),opts.Nobj(2),Nslice);
for m = 1:Nslice
    if z(m) == 0
        Ibf(:,:,m) = sum(lf_stack,3);
    else
        for n = 1:Nangle % BF only
            shift_x = z(m) * tan_theta(n,2);
            shift_y = z(m) * tan_theta(n,1);
            Ibf(:,:,m) = Ibf(:,:,m) + ...
                lf_shift(lf_stack(:,:,n),shift_x,shift_y);
        end
    end
end
Ibf = Ibf/Nangle;

m0 = find(abs(z)==min(abs(z)));
m0 = m0(1);
Imin = min(min(Ibf(:,:,m0)));
Imax = max(max(Ibf(:,:,m0)));


end

