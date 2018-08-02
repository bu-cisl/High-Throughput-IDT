function [imOut] = LPFilter(img,ifgpu)
%%
%LPFilter
%
%Purpose: This is Alex Matlock's copy of the low-pass filtering method used
%in DPC measurements when Dr. Tian was a post-doc. For each image in img
%assuming a 3-D array m x n x p where p = image number, the image is
%filtered by the provided filter h, the background is divided out of the
%image, and the image is rescaled based on the average value of the overall
%image.
%
%Inputs: img - 3-D matrix of images of size m x n x p, p = number of images
%        ifgpu - if GPU is used for computation           
%
%Outputs: - imOut - 3-D stack of filtered images


%Initialize filter, output image stack
if ifgpu
    gimOut = gpuArray(zeros(size(img)));
    imgG = gpuArray(img);

else
    gimOut = zeros(size(img));
    imgG = img;
end

[N1,N2,nImg] = size(img);
h = fspecial('average',[round(N1/3),round(N2/3)]);



%Filter images
for k = 1:nImg
         
    bkgnd = imfilter(imgG(:,:,k),h,'replicate');
    if ifgpu
        gimOut(:,:,k) = imgG(:,:,k) ./ gpuArray(bkgnd);
    else
        gimOut(:,:,k) = imgG(:,:,k) ./ bkgnd;
    end
    
    gimOut(:,:,k) = gimOut(:,:,k)./mean2(gimOut(:,:,k));
    

end %End of image filter loop

if ifgpu
    imOut = gather(gimOut);
else
    imOut = gimOut;
end

end
