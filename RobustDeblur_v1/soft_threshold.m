function ds1=soft_threshold(d, t)
%%%%%%%%%%%%%%%%%%%%%%%%
%Soft thresholding
%
ds=real(d); 
    res = (abs(ds) - t);
    res = (res + abs(res))/2;
    ds_real = sign(ds).*res;
    
ds1 = ds_real;    
% ds=imag(d); 
%    res = (abs(ds) - t);
%    res = (res + abs(res))/2;
%    ds_imag = sign(ds).*res;
    
%    ds1=ds_real + i * ds_imag; 
