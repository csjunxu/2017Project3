function Dec = imFrameDec(A)
    [D,R]=GenerateFrameletFilter(1);
    [co] = size(A,3);
    Dec = cell(1,co);
    for i = 1:co
        Dec{i}=FraDecMultiLevel(A(:,:,i),D,1);
    end
return;
