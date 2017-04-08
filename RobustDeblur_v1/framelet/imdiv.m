function [im]=imdiv(im,h)
    [ny,nx,co] = size(im);
    for i = 1:co
        im(:,:,i) = im(:,:,i)./h;
    end
return;
