function c=isoshrinkage(c, tau)

[co]=length(c);
[D,R]=GenerateFrameletFilter(1);
mu=getwThresh(1,0,1,D);
nD = length(D); 
R=0; % getting weighted local l2 norm R 
for ii = 1:nD-1
    for jj = 1:nD-1
        if ii==1 && jj==1  continue;  end
        R = R + (c{1}{1}{ii,jj}).^2/mu{1}{ii,jj}^2; 
    end
end
R = sqrt(R); 
R = max(R+realmin-tau, 0)./(R+realmin); 
for ii = 1:nD-1
    for jj = 1:nD-1
        if ii==1&&jj==1 continue; end
        c{1}{1}{ii,jj} = c{1}{1}{ii,jj}.*R; 
    end
end
        
        