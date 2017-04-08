function gamma=imFrameCoeffOper(op,alpha,beta)
    [co]=length(alpha);
    [D,R]=GenerateFrameletFilter(1);
    for i = 1:co
        if op=='s' || op=='h'
            mu=getwThresh(beta,0,1,D);
            gamma{i}=CoeffOper(op,alpha{i},mu);
        else
            if isa(beta, 'cell')
                gamma{i}=CoeffOper(op,alpha{i},beta{i});
            else
                gamma{i}=CoeffOper(op,alpha{i},beta);
            end
        end
    end
    if op=='v'
        v = 0;
        for i = 1:co
            v = v+gamma{i};
        end
        gamma = v;
    end
return;
