%Please replace 'model.mat' with a similarly-titled file containing the
%relevant metabolic model. The default model used is that of E. coli
%iAF1260.

load model.mat;
sol=optimizeCbModel(model);
WTgrRate=sol.f;
fs=[];
for i=1:length(model.genes)
    TempModel=deleteModelGenes(model,model.genes(i));
    sol=optimizeCbModel(TempModel);
    fs(end+1)=sol.f;
    disp(['Deleting gene ',num2str(i)])
end

%%
Filtered=[];
for i=1:length(fs)
    if fs(i)>0.1*WTgrRate && fs(i)<0.9*WTgrRate
        Filtered(end+1)=i;
    end
end

%%
serial=num2cell(Filtered');
genes=model.genes(Filtered);
prots=model.proteins(Filtered);
grRatio=fs(Filtered);
grRatio=num2cell(grRatio'/WTgrRate);
showme=[serial,genes,prots,grRatio];
FBAresult=cell2table(showme,'VariableNames',{'Serial Number','Gene Name','Protein Name','grRatio'});
save('FBAresult.mat','FBAresult')
clearvars -except FBAresult