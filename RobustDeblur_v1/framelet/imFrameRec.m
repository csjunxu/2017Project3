function Rec = imFrameRec(C)
    [D,R]=GenerateFrameletFilter(1);
    [co]=length(C);
    for i = 1:co
        Rec(:,:,i)=FraRecMultiLevel(C{i},R,1);
    end
return;
