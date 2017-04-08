function [f, v, w] = APG3_gray(g, BlurOperation, InitialGuess, Transform,  par, option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% robust image deblurring using accelerated proximal method
%
% % Input: 
%       g:             Blurry image
%       BlurOperation: Blur matrix (.mtx)or blur kernel (.psf)
%       InitialGuess:  Initial guess for the clear image (.f)
%       Transform:     Transformation needed in the objective function(.W, .WT, .D, .DT, .F, .FT)
%       par:           Parameters needed, (.lambda1, .lambda2, .lambda3, .L, .beta)
%       option:        Options during the running (.showImg(1 or 0), .silent(1 or 0), .GroundTruthImg)
%
% output:
%   f	-	deblurred result 
%
%
% The optimization model model:
% c*= argmin 0.5*||BW'c + BD'u + F'h-g||_2^2 + 0.5*beta
% ||(I-WW')c||_2^2+||diag(lambda1) c||_1 + ||diag(lambda2) h||_1 +
% ||diag(lambda3) u||_1
%
% B:          Blurring Operator
% W,D,F: Transformation operator
% where W: Wavelet, D: DCT, F: Identity
%
%The minimization problem is solved via the accelerated proximal gradient
%method

%Reference: H. Ji and K. Wang, Robust image deconvolution with an inaccurate blur kernel, IEEE Transactions on Image Processing, 21(4), 1624-1634, Apr. 2012 
%
%Author:  Kang Wang and Hui Ji
%
%Last Revision: 25-May-2014
%



SIZE = [size(g,1), size(g,2)];
if isfield(BlurOperation,'mtx')
    B = @(x)       BlurOperation.mtx * x(:);
    BT = @(x) reshape(x(:)'*BlurOperation.mtx,[SIZE(1), SIZE(2)]);
    PixelsCut = [0 0];
elseif isfield(BlurOperation,'psf')
    PSF  = BlurOperation.psf;
    PSFT = rot90(rot90(PSF));
    PixelsCut = size(PSF);
    g = MirrorExtension(g, SIZE+2*PixelsCut);
    SIZE = [size(g,1), size(g,2)];
    %     B  = @(x) imfilter(x,  PSF,'conv', 'symmetric');
    %     BT = @(x) imfilter( x, PSFT,'conv', 'symmetric');
    % for speed consideration, use the fft2 based image convolution
    PSFext = extendHforConv(PSF, SIZE(1),SIZE(2));
    PSFText = extendHforConv(PSFT, SIZE(1),SIZE(2));
    fPSFext = fft2(PSFext);
    fPSFText = fft2(PSFText);
    B = @(x) ifft2(fft2(x).*fPSFext);
    BT = @(x) ifft2(fft2(x).*fPSFText);
end

Window=[PixelsCut(1)+1, size(g, 1)-PixelsCut(1), PixelsCut(2)+1, size(g, 2)-PixelsCut(2)];


if isfield(InitialGuess, 'f')
    f0 = InitialGuess.f;
    f0 = MirrorExtension(f0, SIZE);
    v0 = zeros(SIZE);
    w0 = v0;
end

if isfield(Transform, 'W')   W  = Transform.W;  end;
if isfield(Transform, 'WT')  WT = Transform.WT; end;
if isfield(Transform, 'F')   F  = Transform.F;  end;
if isfield(Transform, 'FT')  FT = Transform.FT; end;
if isfield(Transform, 'D')   D  = Transform.D;  end;
if isfield(Transform, 'DT')  DT = Transform.DT; end;

if isfield(par,'beta')    beta    = par.beta; end;
if isfield(par,'lambda1') lambda1 = par.lambda1; end;
if isfield(par,'lambda2')
    lambda2 = zeros(size(g));
    lambda2(Window(1):Window(2), Window(3):Window(4)) = par.lambda2;
end

if isfield(par,'lambda3') lambda3 = par.lambda3; end;
if isfield(par,'L')            L  = par.L; end;

if isfield(option, 'nloops')    nloops   = option.nloops; end
if isfield(option, 'showImg')  showImg = option.showImg; end
if isfield(option, 'OutName')  OutName = option.OutName; end

% Initialization
f=f0;
v=zeros(size(f));
c0=W(f);

c1=c0;

h0= F(v0);
h1=h0;

% tmp=D(v);
% u0=zeros(size(tmp));
u0 = D(w0);
u1=u0;

t0=1;
t1=1;

for iter = 1 : nloops
    tic;
    %step 1
    %     dc = c1 + (t0-1)/t1*(c1-c0);
    tmp = imFrameCoeffOper('-', c1, c0);
    tmp = imFrameCoeffOper('*', tmp, (t0-1)/t1);
    
    dc = imFrameCoeffOper('+', c1, tmp);
    dh = h1 + (t0-1)/t1*(h1-h0);
    du = u1 + (t0-1)/t1*(u1-u0);
    
    tmp = imFrameCoeffOper('-', dc, W(WT(dc)));
    tmp = imFrameCoeffOper('*', tmp, beta);
    if isfield(BlurOperation,'psf')
        tmp1 = B(WT(dc) + DT(du)) + FT(dh) -g;
    else
        tmp2 = WT(dc) + DT(du);
        tmp3 = FT(dh)-g;
        tmp1 = A*tmp2(:)+tmp3(:);
    end
    
    
    BTtmp1 = BT(tmp1);
    Grad_c = imFrameCoeffOper('+', W(BTtmp1), tmp);
    Grad_h = F(reshape(tmp1, SIZE)) ; %+ beta * (h1 - F(FT(h1)));
    Grad_u = D(BTtmp1);
    
    
    tmp = imFrameCoeffOper('/', Grad_c, L);
    
    gc = imFrameCoeffOper('-', dc, tmp);
    gh = dh - Grad_h / L;
    
    gu = du - Grad_u / L;
    % step 3
    c0 = c1;
    h0 = h1;
    u0 = u1;
    c1 = imFrameCoeffOper('s', gc, lambda1/L);
    %     c1 = isoshrinkage(gc, lambda1/L);
    h1 = soft_threshold(gh, lambda2/L);
    u1 = soft_threshold(gu, lambda3/L);
    
    
    % step 4
    t0=t1;
    t1=0.5*(1+sqrt(1+4*t0^2));
    
    f = WT(c1);
    v = FT(h1);
    w = DT(u1);
    f(f>1)=1;
    f(f<0)=0;
    c1 = W(f);
    %     f = edgetaper(f, EstPsf);
    %     c1 = W(f);
    
    %     rhs = norm(B(f(:,:,1)+w(:,:,1)) + A(v(:,:,1)) - g(:,:,1), 'fro');
    tend = toc;
    if showImg % show the intermediate results or not
        hold on;
        subplot('position', [0 0 0.8 1]);
        
        imshow(f);
        title(['f, \lambda_1= ', num2str(lambda1)])
        
        subplot('position', [0.8,0.5,0.2,0.5]);
        %         imshow(abs(GroundTruthImg-f)*10);%, 'InitialMagnification', 100);
        imshow(abs(v)*10);
        title(['v, Iteration ', num2str(iter)]);
        
        subplot('position', [0.8,0,0.2,0.5]);
        imshow(abs(w)*10);
        title(['w, \lambda_3=', num2str(lambda3)]);
        hold off;
        drawnow;
        
    end
    
   
    if ~option.silent % print on screeen or not
        if isfield(option, 'GroundTruthImg')
            PSNR = psnr(option.GroundTruthImg, f(Window(1):Window(2), Window(3):Window(4)));
            fprintf('Iter: %3d, time: %.2fs, PSNR:%.3f \n', iter, tend, PSNR);
        else
            fprintf('Iter: %3d, time: %.2fs\n', iter, tend);
        end
    end
end
f = f(Window(1):Window(2), Window(3):Window(4));