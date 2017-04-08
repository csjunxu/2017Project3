function y=psnr(im,denoised_im)

% im = double(im);
% denoised_im = double(denoised_im);
err = im - denoised_im;
%e=norm(err,'fro')^2/prod(size(err));
e = sum(err(:).^2)/numel(err);
y = 20*log10(1/sqrt(e));
