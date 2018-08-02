function [ S ] = SourceFinitePtComp( uled, vled, u, v, dled)
%SourcePtComp calculates the effective source S given
%   Inputs:
%   uled, vled: coordinate of the point source that will be light up in
%               this image
%
%   u,v: spaital frequency axes

S = zeros(size(u));

du = u(1,2) - u(1,1);
dv = v(2,1) - v(1,1);

if dled>=1
    S(((u-uled)/(dled/2*du)).^2+...
        ((v-vled)/(dled/2*dv)).^2<1) = 1;
else
    tmpdist = ((u-uled)/(du/2)).^2+((v-vled)/(dv/2)).^2;
    [~,idx] = min(tmpdist(:));
    S(idx(1)) = 1;
end    

end

