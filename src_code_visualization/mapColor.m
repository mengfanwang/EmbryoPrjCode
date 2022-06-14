function IdxLst=mapColor(connectTable)

colorNum=max(sum(connectTable));
fullSet=1:colorNum;

N=size(connectTable,1);

IdxLst=randi([1 colorNum],1,N);

for iter=1:10
    for i=1:N
        thisCor=IdxLst(i);
        conIdx=connectTable(i,:);
        conCor=IdxLst(conIdx);

        if length(find(conCor==thisCor))>1
            r = setxor(fullSet,unique(conCor));
            if ~isempty(r)
                IdxLst(i)=r(1);
            else
                conCor2=IdxLst(conIdx(conIdx~=i));
                IdxLst(i)=findLeastFreqElement(conCor2);
            end
        end
    end
end

end