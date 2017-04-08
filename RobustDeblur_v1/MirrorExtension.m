function  g=MirrorExtension(g_in, SCALE)
for i=1:size(g_in,3)
    f = g_in(:,:,i);
    x1=ceil((SCALE(1)-size(f, 1))/2);
    x2=SCALE(1)-size(f,1)-x1;
    y1=ceil((SCALE(2)-size(f, 2))/2);
    y2=SCALE(2)-size(f,2)-y1;
    
    g(:,:,i)=[f(x1:-1:1, y1:-1:1) f(x1:-1:1, :) f(x1:-1:1, size(f,2):-1:size(f,2)-y2+1) ;...
        f(:, y1:-1:1) f  f(:, size(f, 2):-1:size(f,2)-y2+1);...
        f(size(f, 1):-1:size(f,1)-x2+1, y1:-1:1) f(size(f, 1):-1:size(f,1)-x2+1,:) f(size(f, 1):-1:size(f,1)-x2+1,size(f, 2):-1:size(f,2)-y2+1)];
end