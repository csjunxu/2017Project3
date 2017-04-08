function demo_details()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The demonstration of robust image deblurring with additional details 
% on change of psnr values during each iteration
%
%Author:  Kang Wang and Hui Ji
%
%Last Revision: 25-May-2014
%

addpath('./Framelet/'); 

% 'cameraman_blurry.png is generated by applying the inear motion kernel
% fspecial('motion', 20,10) on the image 'cameraman.png' with edge chopped
g = im2double(imread('cameraman_blurry.png')); 

% the psf used for deblurring the image 'cameraman_blurry.png' is different
% from the true psf fspecial('motion', 20,10) 
BlurOperation.psf = fspecial('motion', 20,20);

InitialGuess.f = g; 

Transform.W  = @(x) imFrameDec(x); 
Transform.WT = @(x) imFrameRec(x); 
Transform.D  = @(x) dct2(x); 
Transform.DT = @(x) idct2(x); 
Transform.F  = @(x) x; 
Transform.FT = @(x) x; 

par.lambda1 = 0.0003; 
% par.lambda2 is the weight parameter for the outliers. 
% The value of par.lambda2 can be one single value or a matrix with the same size of g
% Single value means uniform constraints for the whole image, while using matrix may you can set different weight for different regions 
% (e.g., you can set the weight for the saturated region to be 0 or very small. The smaller the weight, the more freedom for h.)
% e.g. single value:
% par.lambda2 = 0.0005; 
% e.g., matrix way: 
par.lambda2 = 0.0005*ones(size(g)); 
par.lambda2(g>240/255) = 0; 
par.lambda3 = 0.0003;
par.beta = 1; % suggested value 
par.L = 3; % Larger value will guarantee the convergence while smaller gives faster convergence

option.nloops = 400; %% maximum iteration numbers
option.showImg = 0; 
option.OutName = 'cameraman_out'; 
option.silent = 0; 

% The ground-truth is only served for the purpose of showing the PSNR values during each iteration,
% it is not used in the APG3_gray for deblurring image. 
% One may mark out the line below when used in practice
option.GroundTruthImg = im2double(imread('cameraman.png')); 

[f, ~, ~] = APG3_gray(g, BlurOperation, InitialGuess, Transform,  par, option); 

imwrite(f, [option.OutName, '.png'], 'png'); 