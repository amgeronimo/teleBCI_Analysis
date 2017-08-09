function [] = ShowResultsAllAtOnce()

%This file is the same as ShowResults, except it runs over all the
%participants.  Did this because the cluster wasnt allowing me to run
%multiple parallel jobs



clearvars -except RunNum RNN full_code

runloc = cd;
runloc = runloc(1);



fprintf('*\n***\n*****\n*******\n*********\n***********\n*************\n**************\nFinding P300 Features\n**************\n*************\n***********\n*********\n*******\n*****\n***\n*\n')


[P300] = tBCIFileNames();
Names = P300.Names; Session = P300.Session; Run = P300.Run;
CapVersion = P300.CapVersion;
%Feedback = P300.Feedback; FreeSpell = P300.FreeSpell; Type = P300.Type;







for RunNum = 1:length(Names)
    
    
    
    %% Impedance Extractor (P300)
    tempName = Names{RunNum};
    tempSession = Session{RunNum};
    
    
    
    clear IMP
    ElecNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};
    kk = 1;
    for sess = 1:length(tempSession)
        disp(['IMP_' tempName '_S' tempSession{sess}])
        
        if runloc == 'C'
            d = dir(['data\' tempName tempSession{sess} '\impedancefile*']);
            if isempty(d)
                fid = -1;
            else
                fid = fopen(['data\' tempName tempSession{sess} '\' d(end).name],'r');
                if fid == -1
                    fid = fopen(['P:\ALS Proj Data\' tempName ...
                        tempSession{sess} '\' d.name],'r');
                end
            end
            
        else
            d = dir(['/gpfs/work/a/amg5106/ALS_P300/Data/' tempName ...
                tempSession{sess} '/impedancefile*']);
            fid = fopen(['/gpfs/work/a/amg5106/ALS_P300/Data/' tempName ...
                tempSession{sess} '/' d.name],'r');
        end
        
        if fid == -1
            disp('Impedance File not found');
        else
            itt=textscan(fid,'%s','delimiter','\n');
            fclose(fid)
            InputText = itt{1};
            clear e_ind e_loc
            for jj = 1:length(ElecNames)
                for ll = 1:length(InputText)
                    e_ind(jj,ll) = ~isempty(strfind(InputText{ll},ElecNames{jj}));
                end
                e_loc{jj} = find(e_ind(jj,:)==1);
            end
            e_loc_end = cellfun(@(x) x(end),e_loc);
            IMPS = InputText(e_loc_end);
            expr = [':' '(.*?)' 'kOhm'];
            tok = cellfun(@(x) regexp(x,expr,'tokens'),IMPS);
            IMP(kk,:)=cellfun(@(x) str2num(cell2mat(x)),tok,'UniformOutput',false);
            labelP{kk} = [tempName tempSession{sess}];
            
            %         [CData, NData, Spec] =  ShowMeTheP300(tempName,tempSession{sess},...
            %             tempRun(ii),0,2,0, full_code);
            %
            %         CD(kk,:,:) = mean(CData,3);
            %         ND(kk,:,:) = mean(NData,3);
            %         SD(kk,:,:) = Spec;
            kk = kk+1;
        end
    end
    IMP = cell2mat(IMP);
    save(['results/' tempName '_P300_TrainerImpedance.mat'],'IMP','labelP');
    
    
    %% ShowMeTheP300
    
    %To generate Name_P300 feats files, have Centerfocus==0, FreeSpell==0,
    %Speed3x == 0, and change the name accordingly at the end
    %To generate the Name_CF_P300... files, change Centerfocus==1 and the
    %savename
    %To generate the Name_fast_P300... files, change Speed3x==0 and the
    %savename
    
    Rerefs = [0 1];
    ARTs = [0 1];
    
    
    
    
    
    
    
    %
    %     Indices = FreeSpell == 0;
    %     Indlocs = find(Indices==1);
    %     tempNames = Names(Indices);
    %     tempSession = Session(Indices);
    %     tempRun = Run(Indices);
    %     tempFeedback = Feedback(Indices);
    %     nms = unique(tempNames);
    %
    %
    %
    %     tempRuns = Run{RunNum}{ss}
    %     tempFree = cell2mat(FreeSpell{RunNum}{ss});
    %     tempFeedback = cell2mat(Feedback{RunNum}{ss});
    %     tempType = cell2mat(Type{RunNum}{ss});
    %
    %
    %     Indices = tempFree == 0 & tempFeedback == 0;
    %     Indlocs = find(Indices==1);
    %             tempNames = Names(Indices);
    %             tempSession = Session(Indices);
    %             tempRun = Run(Indices);
    %             tempFeedback = Feedback(Indices);
    %
    %             nms = unique(tempNames);
    %             runstodo = find(strcmp(tempNames,nms(RunNum)));

    for ART = 1:length(ARTs)
    for RV = 1:length(Rerefs)
        
        %                 subj_files = strcmp(tempNames,nms{subj})==1;
        %                 ses = unique(tempSession(subj_files));
        clear C_m N_m SPEC_m C_s N_s SPEC_s C_e_m C_l_m...
            N_e_m N_l_m C_e_s C_l_s N_e_s N_l_s Ntime Ctime...
            C_keep N_keep
        C_FULL = []; N_FULL = []; Sess_FULL = [];
        C_FULLe = []; N_FULLe = []; C_FULLl = []; N_FULLl = [];
        Csestrl_FULL = []; Nsestrl_FULL = [];
        AIc_ALL = []; AIn_ALL = [];
        TW_ALL= [];
        GAZE_FULL = [];
        FaceData_FULL = cell(1,24);
        
        for sess = 1:length(tempSession)
            tempRuns = Run{RunNum}{sess};
            
            thesefiles = 1:length(tempRuns);
            
            C_ALL = [];
            N_ALL = [];
            SPEC_ALL = zeros(251,8,length(thesefiles));
            BSPEC_ALL = zeros(251,8,length(thesefiles));
            CEOG_ALL = [];
            NEOG_ALL = [];
            C_EARLY = []; C_LATE = [];
            N_EARLY = []; N_LATE = [];
            
            for ii = 1:length(thesefiles)
                %      if strcmp(Names{ii},'P11')==1
                disp([tempName ' ' num2str(RV) ' ' num2str(ART) ' ' tempSession{sess} ' ' tempRuns{thesefiles(ii)}])
                close all
                [CData, NData, ~, ~, Nsq, AIc, AIn, TargWait, GAZE, FaceData] = ...
                    ShowMeTheP300_v2(tempName,tempSession(sess),...
                    tempRuns(thesefiles(ii)),Rerefs(RV),ARTs(ART),0,CapVersion(RunNum));
                
%                 if sum(sum(sum(isnan(NData))))~=0
%                     hey = 1
%                 end
                if  Rerefs(RV)==0
                    [Spec, BSpec, ff] = ShowMeTheSpectra(tempName,tempSession{sess},...
                        tempRuns(thesefiles(ii)),Rerefs(RV),0,0);
                    if ~isempty(Spec)
                        SPEC_ALL = cat(3,SPEC_ALL,Spec);
                    else
                        SPEC_ALL = cat(3,SPEC_ALL,NaN(251,8));
                    end
                    if ~isempty(BSpec)
                        BSPEC_ALL = cat(3,BSPEC_ALL,BSpec);
                    else
                        BSPEC_ALL = cat(3,BSPEC_ALL,NaN(251,8));
                    end
                end
                
                
                GAZE_FULL{sess,ii} = GAZE;
                Csestrl_FULL = [Csestrl_FULL;...
                    [repmat(str2num(tempSession{sess}),size(CData,3),1)...
                    repmat(ii,size(CData,3),1) zeros(size(CData,3),1)]];
                Nsestrl_FULL = [Nsestrl_FULL;...
                    [repmat(str2num(tempSession{sess}),size(NData,3),1)...
                    repmat(ii,size(NData,3),1) zeros(size(NData,3),1)]];
                C_FULL = cat(3,C_FULL,CData);
                N_FULL = cat(3,N_FULL,NData);
                
                
                
                C_ALL = cat(3,C_ALL,nanmean(CData,3));
                N_ALL = cat(3,N_ALL,nanmean(NData,3));
                
                %                 if Rerefs(RV)==0
                %                     for ff = 1:24
                %                         if ~isempty(FaceData.C{ff})
                %                             FaceData_FULL{ff} = cat(3,FaceData_FULL{ff},nanmean(FaceData.C{ff},3));
                %                         end
                %                     end
                %                 end
                
                
                AIc_ALL = [AIc_ALL AIc];
                AIn_ALL = [AIn_ALL AIn];
                
                TW_ALL = [TW_ALL TargWait];
                
                
                
                
                
            end
            C_keep = size(C_FULL,3);
            N_keep = size(N_FULL,3);
            
            
            
            if  Rerefs(RV)==0
                
                if length(size(SPEC_ALL))==3 %More than one run
                    SPEC_m(sess,:,:) = single(nanmean(SPEC_ALL,3));
                    SPEC_s(sess,:,:) = single(nanstd(SPEC_ALL,[],3));
                end
                if  length(size(BSPEC_ALL))==3
                    BSPEC_m(sess,:,:) = single(nanmean(BSPEC_ALL,3));
                    BSPEC_s(sess,:,:) = single(nanstd(BSPEC_ALL,[],3));
                end
                if length(size(SPEC_ALL))==2 %Just one run
                    if ~isempty(SPEC_ALL)
                        SPEC_m(sess,:,:) = single(SPEC_ALL);
                    else
                        SPEC_m(sess,:,:) = NaN(size(SPEC_ALL));
                    end
                    SPEC_s(sess,:,:) = NaN(size(SPEC_ALL));
                end
                if  length(size(BSPEC_ALL))==2
                    if ~isempty(BSPEC_ALL)
                        BSPEC_m(sess,:,:) = single(BSPEC_ALL);
                    else
                        BSPEC_m(sess,:,:) = NaN(size(BSPEC_ALL));
                    end
                    BSPEC_s(sess,:,:) = NaN(size(BSPEC_ALL));
                end
            else
                SPEC_m = []; SPEC_s = []; BSPEC_m = []; BSPEC_s = [];
            end
            
            
            if length(size(C_ALL))==3 %more than one run
                C_m(sess,:,:) = single(nanmean(C_ALL,3));
                N_m(sess,:,:) = single(nanmean(N_ALL,3));
                C_s(sess,:,:) = single(nanstd(C_ALL,[],3));
                N_s(sess,:,:) = single(nanstd(N_ALL,[],3));
                
                
            elseif length(size(C_ALL))==2
                if ~isempty(thesefiles) %There was one run
                    C_m(sess,:,:) = single(C_ALL);
                    N_m(sess,:,:) = single(N_ALL);
                    C_s(sess,:,:) = single(C_ALL);
                    N_s(sess,:,:) = single(N_ALL);
                    
                    
                    
                else %no runs of this type
                    sizeofdata = [188 9];  %This is not good, having
                    %predefine the size of the empty matrices, but w/e.
                    sizeofspec = [251 9];
                    C_m(sess,:,:) = single(NaN(sizeofdata));
                    N_m(sess,:,:) = single(NaN(sizeofdata));
                    %                             SPEC_m(sess,:,:) = single(NaN(sizeofspec));
                    %                             BSPEC_m(sess,:,:) = single(NaN(sizeofspec));
                    C_s(sess,:,:) = single(NaN(sizeofdata));
                    N_s(sess,:,:) = single(NaN(sizeofdata));
                    %                             SPEC_s(sess,:,:) = single(NaN(sizeofspec));
                    %                             BSPEC_s(sess,:,:) = single(NaN(sizeofspec));
                    
                end
                
            end
        end
        size(C_FULL)
        cfm = nanmean(C_FULL,3);
        cfs = nanstd(C_FULL,[],3);
        nfm = nanmean(N_FULL,3);
        nfs = nanstd(N_FULL,[],3);
        
        
        
        AIC = sum(AIc_ALL,2);
        AIN = sum(AIn_ALL,2);
        
%         %Try doing a weighted average using the time
%         %preceding the target.  Weight those trials with a
%         %long wait time as higher
%         TWkeep = TW_ALL<4*256 & TW_ALL>.25*256; %Keep trials
%         %which were preceded by another flash no less than
%         %.25 seconds and no more than 4 seconds.
%         TW_ALL = TW_ALL(TWkeep);
%         
%         %                         %For image figure in paper P300_W_hist.pdf
%         %                          [bhist bloc] = hist(TW_ALL,14);
%         %                          bhist = [bhist;bhist];
%         %                          bhist(1,6:8) = 0;
%         %                          figure('renderer','painters');
%         %                          hb = bar3(bloc/256,bhist'); view([-70 30]);
%         %                          zlabel('Trial Count'); ylabel('Inter-target interval (seconds)');
%         %                          xlim([0 3]);ylim([0 4]);
%         %                          set(gca,'Xticklabel',[]);
%         %                          set(hb(2),'facecolor',[1 0 0]);
%         %                          set(hb(1),'facecolor',[0 0 1]);
%         
%         
%         C_FULLW = C_FULL(:,:,TWkeep);
%         
%         %Then remove flashes around the median time
%         TWkeep2 = TW_ALL ~= 448 & TW_ALL ~=512 & TW_ALL ~= 576;
%         TW_ALL = TW_ALL(TWkeep2);
%         C_FULLW = C_FULLW(:,:,TWkeep2);
%         
%         %Was doing this, but it ended up just reducing the
%         %overall VEP all over.
%         %                         %Weight trials so that they are worth double
%         %                         %for short flash time and half for median flash
%         %                         %times
%         %                         for ll = 128:64:448
%         %                             C_FULLW(:,:,TW_ALL==ll) = C_FULLW(:,:,TW_ALL==ll)*(512/ll)*.5;
%         %                         end
%         %                         for ll = 512:64:960
%         %                             C_FULLW(:,:,TW_ALL==ll) = C_FULLW(:,:,TW_ALL==ll)*(ll/512)*.5;
%         %                         end
%         cfmW = mean(C_FULLW,3);
%         cfsW = std(C_FULLW,[],3);
%         CW_keep = size(C_FULLW,3);
%         clear C_FULLW;
        
        %                 %Average ans std the first, second, third runs over
        %                 %sessions
        %                 for sess = 1:max(Nsestrl_FULL(:,2))
        %                     Ntime.dm(sess,:,:) = mean(N_FULL(:,:,Nsestrl_FULL(:,2)==sess),3);
        %                     Ntime.ds(sess,:,:) = std(N_FULL(:,:,Nsestrl_FULL(:,2)==sess),[],3);
        %                     Ntime.count(sess) = sum(Nsestrl_FULL(:,2)==sess);
        %                     Ctime.dm(sess,:,:) = mean(C_FULL(:,:,Csestrl_FULL(:,2)==sess),3);
        %                     Ctime.ds(sess,:,:) = std(C_FULL(:,:,Csestrl_FULL(:,2)==sess),[],3);
        %                     Ctime.count(sess) = sum(Csestrl_FULL(:,2)==sess);
        %                 end
        
        
        
     save(['results/' tempName '_P300feats_Ref' num2str(Rerefs(RV)) '.mat'],'C_m','N_m',...
        'C_s','N_s',...
        'C_keep','N_keep','cfm','cfs',...
        'nfm','nfs','AIC','AIN',...
        'GAZE_FULL','FaceData_FULL',...
        'SPEC_m','SPEC_s','BSPEC_m','BSPEC_s','ff');    
    end
   
    %Took out, but can put back in selectively for ART = 0 and ART = 2
    %cfem','cfes','nfem','nfes','cflm',...
    %    'cfls','nflm','nfls'
    %'SPEC_m','SPEC_s','N_e_m','N_e_s','N_l_m','N_l_s','C_e_m','C_e_s','C_l_m','C_l_s'
    
    
    
    end   
    
    


end
