changed_M=Media(:,1);
changed_M=strcat('R_',regexprep(changed_M,'\(u)','_e'));
j=0;
Pos=length(changed_M);
% i = 1: length(Reference)
for i = 1: length(Reference)
    Temp=Reference{i}.model;
    Temp.lb(findExcRxns(Temp))=0;
    Temp=changeRxnBounds(Temp,changed_M,-10,'l');
    [a,b]=findMinimalMedia(Temp);
    if ~isempty(b)
        changed_M=[changed_M;b]
    end
    optimizeCbModel(Temp)
end