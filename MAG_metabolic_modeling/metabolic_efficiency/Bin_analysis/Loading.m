
%% Adding the growth media
% fname = 'Media_Ready.json';
% Complete =jsondecode(fileread(fname));
% 
% 
% fname = 'Minerals.json';
% Minerals = jsondecode(fileread(fname));
% Temp_CM=fieldnames(Minerals.CM);
% for i=1:length(Temp_CM)
% Mineral{i,1}=Minerals.CM.(Temp_CM{i})
% end
% Mineral(:)=regexprep(Mineral(:),'-','_');
% Mineral(:)=strcat('EX_',regexprep(Mineral(:),'\[u]','(u)'));
% Minerals=Mineral;
% clear Mineral
% % 
% % clear raw str fname fid ans
% % fname = 'Vitamins.json';
% % fid = fopen(fname);
% % raw = fread(fid);
% % str = char(raw');
% % fclose(fid);
% % Vitamins = JSON.parse(str);
% %     
% % clear raw str fname fid ans
% fname = 'Lipid.json';
% Lipid = jsondecode(fileread(fname));
% Temp_CM=fieldnames(Lipid.CM);
% for i=1:length(Temp_CM)
% Lipids{i,1}=Lipid.CM.(Temp_CM{i})
% end
% Lipids(:)=regexprep(Lipids(:),'-','_');
% Lipids(:)=strcat('EX_',regexprep(Lipids(:),'\[u]','(u)'));
% % fname = 'Fiber.json';
% % fid = fopen(fname);
% % raw = fread(fid);
% % str = char(raw');
% % fclose(fid);
% % Fiber = JSON.parse(str);
% % 
% % clear raw str fname fid ans
% 
% fname = 'Fiber.json';
% Fibers = jsondecode(fileread(fname));
% Temp_CM=fieldnames(Fibers.CM);
% for i=1:length(Temp_CM)
% Fiber{i,1}=Fibers.CM.(Temp_CM{i})
% end
% Fibers=Fiber;
% Fibers(:)=regexprep(Fibers(:),'-','_');
% Fibers(:)=strcat('EX_',regexprep(Fibers(:),'\[u]','(u)'));
% 
% clear Fiber
% 
% 
% 
% 
% 
% 
% 
% 
% 
% fname = 'Carbs.json';
% Carbs = jsondecode(fileread(fname));
% Temp_CM=fieldnames(Carbs.CM);
% for i=1:length(Temp_CM)
% Carb{i,1}=Carbs.CM.(Temp_CM{i})
% end
% Carbs=Carb;
% 
% Carbs(:)=regexprep(Carbs(:),'-','_');
% Carbs(:)=strcat('EX_',regexprep(Carbs(:),'\[u]','(u)'));
% 
% Carbs_D=regexprep(Carbs,'__D','__L');
% Carbs=[Carbs;Carbs_D];
% Carbs=unique(Carbs);
% 
% clear Carb Carbs_D
% % clear raw str fname fid ans
% % 
% % fname = 'PolyAmine';
% % fid = fopen(fname);
% % raw = fread(fid);
% % str = char(raw');
% % fclose(fid);
% % PolyAmine = JSON.parse(str);
% % clear raw str fname fid ans
% % fname = 'Water.json';
% % fid = fopen(fname);
% % raw = fread(fid);
% % str = char(raw');
% % fclose(fid);
% % Water = JSON.parse(str);
% % clear raw str fname fid ans
% fname = 'AAs.json';
% AAs = jsondecode(fileread(fname));
% Temp_CM=fieldnames(AAs.CM);
% for i=1:length(Temp_CM)
% AA{i,1}=AAs.CM.(Temp_CM{i})
% end
% AAs=AA;
% 
% 
% AAs(:)=regexprep(AAs(:),'-','_');
% AAs(:)=strcat('EX_',regexprep(AAs(:),'\[u]','(u)'));
% AAs_D=regexprep(AAs,'__L','__D');
% AAs=[AAs;AAs_D]
% clear AA AAs_D
% 
% 
% 
% 
% fname = 'Vitamins.json';
% Vitas = jsondecode(fileread(fname));
% Temp_CM=fieldnames(Vitas.CM);
% for i=1:length(Temp_CM)
% Vitamins{i,1}=Vitas.CM.(Temp_CM{i})
% end
% 
% Vitamins(:)=regexprep(Vitamins(:),'-','_');
% Vitamins(:)=strcat('EX_',regexprep(Vitamins(:),'\[u]','(u)'));
% Vitamins_D=regexprep(Vitamins,'__L','__D');
% Vitamins=unique([Vitamins;Vitamins_D])
% clear Vitas Vitamins_D
% 

% clear raw str fname fid ans
% fname = 'Minerals.json';
% fid = fopen(fname);
% raw = fread(fid);
% str = char(raw');
% fclose(fid);
% Minerals = JSON.parse(str);
% clear raw str fname fid ans
% Common={};
% Temp_CM=fieldnames(Minerals.CM);
% for i=1:length(fieldnames(Minerals.CM))
% Common{i,1}=Minerals.CM.(Temp_CM{i})
% Common{i,2}=Minerals.MW.(Temp_CM{i})
% end
% Pos=length(fieldnames(Minerals.CM));
% Temp_CM=fieldnames(Vitamins.CM);
% for i=1:length(fieldnames(Vitamins.CM))
% Common{Pos+i,1}=Vitamins.CM.(Temp_CM{i})
% Common{Pos+i,2}=Vitamins.MW.(Temp_CM{i})
% end
% Pos=size(Common,1);
% Common{Pos+1,1}=Water.CM;
% Common{Pos+1,2}=Water.MW;

% Pos=size(Common,1);
% Lipid_Rich_Medium=Common;
% Temp_CM=fieldnames(Lipid.CM);
% for i=1:length(fieldnames(Lipid.CM))
%     Lipid_Rich_Medium{Pos+i,1}=Lipid.CM.(Temp_CM{i});
%     Lipid_Rich_Medium{Pos+i,2}=Lipid.MW.(Temp_CM{i});
% end
% 
% AA_Rich_Medium=Common;
% Temp_CM=fieldnames(AAs.CM);
% for i=1:length(fieldnames(Lipid.CM))
%     AA_Rich_Medium{Pos+i,1}=AAs.CM.(Temp_CM{i});
%     AA_Rich_Medium{Pos+i,2}=AAs.MW.(Temp_CM{i});
% end
% 
% Carb_Rich_Medium=Common;
% Temp_CM=fieldnames(Carbs.CM);
% for i=1:length(fieldnames(Carbs.CM))
%     Carb_Rich_Medium{Pos+i,1}=Carbs.CM.(Temp_CM{i});
%     Carb_Rich_Medium{Pos+i,2}=Carbs.MW.(Temp_CM{i});
% end
% 
% Fiber_Rich_Medium=Common;
% Temp_CM=fieldnames(Fiber.CM);
% for i=1:length(fieldnames(Fiber.CM))
%     Fiber_Rich_Medium{Pos+i,1}=Fiber.CM.(Temp_CM{i});
%     Fiber_Rich_Medium{Pos+i,2}=Fiber.MW.(Temp_CM{i});
% end

% Medium.Lipid=Lipid_Rich_Medium;
% Medium.AAs=AA_Rich_Medium;
% Medium.Carbs=Carb_Rich_Medium;
% Medium.Fiber=Fiber_Rich_Medium;

% clearvars -except Medium Common
% fname = 'Complete_Media.json';
% Complete_Media =jsondecode(fileread(fname));
% Fields=fieldnames(Complete_Media.CM);
% for i=1:length(Fields)
%     Media{i,1}=Complete_Media.CM.(Fields{i})
%     if ischar(Complete_Media.GramCpermole.(Fields{i}))
%         Media{i,2}=str2num(Complete_Media.GramCpermole.(Fields{i}))*12;
%     else
%         Media{i,2}=0
%     end
% end

% Pos=length(Media(:,1));
% Media(Pos+1:Pos+7,1)={'n2[u]',
% 'n2o[u]',
% 'nh4[u]',
% 'no[u]',
% 'no2[u]',
% 'no3[u]',
% 'o2[u]'}; 
% 
% Media(Pos+1:Pos+7,2)={0;
% 0;
% 0;
% 0;
% % 0;
% % 0;0}; 
% Media(:,1)=regexprep(Media(:,1),'-','_');
% Media(:,1)=strcat('EX_',regexprep(Media(:,1),'\[u]','(u)'));

load Media 

AAs=find(ismember(Media(:,2),{'AA'       }));
Carbs=find(ismember(Media(:,2),{'Carb'       }));
Lipid=find(ismember(Media(:,2),{'Lipid'       }));

General=[AAs;Carbs;Lipid];
All=1:size(Media,1);
All=All';
General=find(~ismember(All,General));

% Media=load('Complete_Media_Corrected.mat').changed_M;
% Media=regexprep(Media,'_e','(u)');
% Media=regexprep(Media,'R_','');
% Media(:,2)=num2cell(ones(length(Media),1));


fname = 'Candidates.json';
data =jsondecode(fileread(fname));
FNs=fieldnames(data)

fname = 'Rel_Abunds.json';
Relative_Abundances = jsondecode(fileread(fname));

%% Loading the Models
cd './carveme_xml_files/'
Set=fieldnames(data);
UNION={}
for i=1:length(Set);
UNION=union(UNION,data.(Set{i}));
end
Reference={};
for i=1:length(UNION)
    Temp=load(UNION{i});
    Reference{i}=Temp;
    i
end
cd ..


%% Building the community model 
for i=300:length(FNs)
    Models=data.(FNs{i});
    Rel_Abunds=Relative_Abundances.(FNs{i});
    for j=1:length(Models)
        Abund_Model{j}=Reference{find(ismember(UNION,Models{j}))}.model;
        j
        Abund_Model{j}.mets = regexprep(Abund_Model{j}.mets, 'C_', '');
        nameTagsModel{j} = Abund_Model{j}.modelID;
    end
    disp('/nAll Models loaded successfully!')
    options.spBm=cell(1,length(Abund_Model));
    options.spBm(:)={'Growth'};
    options.spATPM=cell(1,length(Abund_Model));
    options.spATPM(:)={'R_ATPM'};
    options.sepUtEx = false;
    MAGCom = createCommModel(Abund_Model, options);
    TempCom=MAGCom;
    TempCom.lb(find(findExcRxns(TempCom)))=0;
    Biomass_Array=MAGCom.rxns(find(MAGCom.c));
    for T=2:length(Biomass_Array)
        RXNCELL{1}=Biomass_Array{T};
        RXNCELL{2}=Biomass_Array{1};
       TempCom =addRatioReaction(TempCom, RXNCELL, [Rel_Abunds(T),Rel_Abunds(1)]);
    end
    
    TempCom=changeRxnBounds(TempCom,Media(General,1),-10000,'l');
    Secondary=zeros(length(TempCom.rxns),1);
    for k=1:length(Media(:,1))
          ind=find(ismember(TempCom.rxns,Media(k,1)));
          if ind 
                 Secondary(ind)=1;
          end
    end
    Secondary(find(TempCom.c))=1;
    Carbs_TempCom=TempCom;
    Carbs_TempCom=changeRxnBounds(Carbs_TempCom,Media(Carbs,1),-100,'l');
    Carbs_TempCom=changeRxnBounds(Carbs_TempCom,Media(Lipid,1),-0.01,'l');
    Carbs_TempCom=changeRxnBounds(Carbs_TempCom,Media(AAs,1),-0.01,'l');
    Sol=optimizeCbModel(Carbs_TempCom,'max')
    Uptake_info(i,1).values= Sol.x(find(Secondary~=0));
    Uptake_info(i,1).rxns= TempCom.rxns(find(Secondary~=0));
    
    
    Lipids_TempCom=TempCom;
    Lipids_TempCom=changeRxnBounds(Lipids_TempCom,Media(Lipid,1),-100,'l');
    Lipids_TempCom=changeRxnBounds(Lipids_TempCom,Media(Carbs,1),-0.01,'l');
    Lipids_TempCom=changeRxnBounds(Lipids_TempCom,Media(AAs,1),-0.01,'l');
    Sol=optimizeCbModel(Lipids_TempCom,'max')
    Uptake_info(i,2).values= Sol.x(find(Secondary~=0));
    Uptake_info(i,2).rxns= TempCom.rxns(find(Secondary~=0));
   
    AAs_TempCom=TempCom;
    AAs_TempCom=changeRxnBounds(AAs_TempCom,Media(AAs,1),-100,'l');
    AAs_TempCom=changeRxnBounds(AAs_TempCom,Media(Lipid,1),-0.01,'l');
    AAs_TempCom=changeRxnBounds(AAs_TempCom,Media(Carbs,1),-0.01,'l');
    Sol=optimizeCbModel(AAs_TempCom,'max')
    
%     printUptakeBoundCom(Carbs_TempCom)
%     printUptakeBoundCom(Lipids_TempCom)
%     printUptakeBoundCom(AAs_TempCom)
    Uptake_info(i,3).values= Sol.x(find(Secondary~=0));
    Uptake_info(i,3).rxns= TempCom.rxns(find(Secondary~=0));


        fprintf('\nSample %d finished successfully!',i)
end
            
    
save Uptake_info
    
    
    
cd ..

