%%




%teleBCI Analysis
close all
clear all

cdp = cd;
[~, sl] = regexp(cdp,'\Users\');
se = regexp(cdp,'Box Sync');
usrnm = cdp(sl+1:se-2);


%pathToFile = ['C:\Users\' usrnm '\Box Sync\Hershey_2016\ALSA teleBCI\Analysis\'];
% %Load data from REDCap.
[PItems, PData] = ...
    ReadFromExcelFile_teleBCI('TeleBCI_DATA_2017-08-07_0951');

IND.code = 1;
IND.consent = 70;

P.fullcode = PData(:,IND.code);
tmp = PData(:,IND.consent);
cons = cellfun(@(x) sum(~isnan(x))>0, tmp);
P.code = unique(P.fullcode);
tmp = ismember(P.code,P.fullcode(cons));
P.code = P.code(tmp);
for i = 1:length(P.code)
    tmp = find(strcmp(PData(:,1),P.code{i}));
    tmp2 = find(~cellfun(@isempty,regexp(PData(tmp,2),'initial')));
    P.initial{i} = tmp(tmp2);
    tmp2 = find(~cellfun(@isempty,regexp(PData(tmp,2),'telebci')));
    P.tbcisess{i} = tmp(tmp2);
    tmp2 = find(~cellfun(@isempty,regexp(PData(tmp,2),'clinical')));
    P.interaction{i} = tmp(tmp2);
    tmp2 = find(~cellfun(@isempty,regexp(PData(tmp,2),'final')));
    P.final{i} = tmp(tmp2);
end
% P.sessnum = [diff(P.firstsess); (length(P.fullcode)+1)-P.firstsess(end)];
% P.firstsess = P.firstsess(tmp);
% P.sessnum = P.sessnum(tmp);
% P.tbcisess = P.sessnum-1;
% %Remove sessions that were the interaction
% tmp2 = regexp(PData(P.firstsess+P.tbcisess+1,2),'clinical');
% P.interaction = NaN(length(P.code),1);
% for i= 1:length(P.code)
%     if ~isempty(tmp2{i})
% P.interaction(i) = P.firstsess(i)+P.tbcisess(i)+1;
%     end
% end
% tmp3 = ~cellfun(@isempty,tmp2);
% P.tbcisess(tmp3) = P.tbcisess(tmp3)-1;
% %Remove sessions that were the final study session
% tmp2 = regexp(PData(P.firstsess+P.tbcisess,2),'final');
% P.finalsess = NaN(length(P.code),1);
% for i= 1:length(P.code)
%     if ~isempty(tmp2{i})
% P.finalsess(i) = P.firstsess(i)+P.tbcisess(i);
%     end
% end
% tmp3 = ~cellfun(@isempty,tmp2);
% P.tbcisess(tmp3) = P.tbcisess(tmp3)-1;


IND.DOB = 72;
IND.sex = 73;
IND.ed = 74;
IND.marital = 75;
IND.DOSO = 76;

IND.date_first = 26;
IND.alsfrs = 81:94;
IND.alsfrs_fu = 934:947;
IND.ecas = 96:110;
IND.alssqol = 113:119;
IND.alssqol_fu = 950:956;
IND.hosu = 120:635;
IND.hosu_fu = 957:1472;
IND.adtpa = 636:719;
IND.adtpa_fu = 1473:1556;
IND.vision = 720:738;

IND.log_date = 741;
IND.log_comm = 743:745;
IND.log_help = 746:757;
IND.log_num = 758;
IND.log_rate = 759:764;
IND.log_speller = 765:768;
IND.log_access = 769:771;
IND.log_conf = 772:781;
IND.log_other = 782:785;

IND.rlog_date = 787;
IND.rlog_comm = 789:791;
IND.rlog_help = 792:803;
IND.rlog_time = 804:813;
IND.rlog_conf = 814:823;
IND.rlog_speller = 824:827;
IND.rlog_access = 828:830;

IND.intx_id = 839;
IND.intx_date = 840;
IND.intx_assess = 845:910;
IND.intx_nassess = 913:930;


%structfun(@(x) x+18,IND,'UniformOutput',false)

keepdata = cellfun(@(x) sum(~isnan(x))>1, (PData(cell2mat(P.initial),IND.rlog_date)));
P.code = P.code(keepdata);
P.initial = P.initial(keepdata);
P.tbcisess = P.tbcisess(keepdata);
P.interaction = P.interaction(keepdata);
P.final = P.final(keepdata);

P.interaction(cellfun(@isempty,P.interaction))= {NaN};
P.final(cellfun(@isempty,P.final))= {NaN};

%% Demographics

fstemp = cell2mat(P.initial);
P.DOB = datestr(PData(fstemp,IND.DOB),'mm/dd/yyyy');
% P.DOSO = datestr(PData(P.firstsess,IND.DOSO),'mm/dd/yyyy');
P.sex = cell2mat(PData(fstemp,IND.sex));
P.ed = cell2mat(PData(fstemp,IND.ed));
P.marital = cell2mat(PData(fstemp,IND.marital));
P.DOS1 = datestr(PData(fstemp,IND.rlog_date),'mm/dd/yyyy');
P.age = floor(daysact(cat(1,P.DOB),cat(1,P.DOS1))/365.25);
P.alsfrs = cell2mat(PData(fstemp,IND.alsfrs(end)));
P.alsfrs_fu = cell2mat(PData(fstemp+cellfun(@length,P.tbcisess)+...
    ~cellfun(@isnan,P.interaction)+~cellfun(@isnan,P.final),IND.alsfrs_fu(end)));
P.alsfrs_fu(isnan(P.alsfrs_fu)) = P.alsfrs(isnan(P.alsfrs_fu));
P.vision = cell2mat(PData(fstemp,IND.vision))

GenCode = {'F','M'};
EdCode = {'No High School','Some High School','High School/GED','Voc. or Tech School',...
    'Some College','College Degree','Graduate Degree'};

SpellerTYPE = {'V','V/A','V','A','V','V'};  %This needs to be manually changed 


%Need to add DOSO when dates are available
DemT = table(P.age,GenCode(P.sex)',EdCode(P.ed)',P.alsfrs, SpellerTYPE');
DemT.Properties.RowNames = P.code;
DemT.Properties.VariableNames = {'Age','Gender','Education','ALSFRSR','SpellerType'}


%% Researcher and Participant confidence.
ActivityLabels = {'Computer','AdobeC','Eye','Cap','Elec','BCI2000','Imp','P300','Update','Clean'};
hh = figure('Position',[100 100 1400 300]);
cmp = jet(10);
clear researcherconf participantconf
for i = 1:length(P.code)
researcherconf(:,:,i) = [-cell2mat(PData(P.initial{i} +[0:length(P.tbcisess{i})],IND.rlog_conf));
    NaN(8-length(P.tbcisess{i}),length(IND.rlog_conf))];
researcherconf(:,:,i) = researcherconf(:,:,i) + repmat(linspace(-.02,.02,10),size(researcherconf,1),1)

participantconf(:,:,i) = [-cell2mat(PData(P.initial{i} +[0:length(P.tbcisess{i})],IND.log_conf));
     NaN(8-length(P.tbcisess{i}),length(IND.log_conf))];
participantconf(:,:,i) = participantconf(:,:,i) + repmat(linspace(-.02,.02,10),size(participantconf,1),1)
end
SUPERPLOT(1,4,1,2,.2,.1,.05);
for i = 1:10
    plot(nanmean(researcherconf(:,i,:),3),'LineWidth',1,'Color',cmp(i,:)); hold on;
end
ylim([-4.1 -.7]); xlim([1 size(researcherconf,1)]); title('Researcher');
plot(nanmean(nanmean(researcherconf,3),2),'LineWidth',3,'Color','k');
set(gca,'Yticklabel',[],'FontSize',12)
textbank =  {{'Not at all','Confident'},{'A Little','Confident'},...
    {'Somewhat','Confident'},{'Very','Confident'}};
for i = -4:1:-1
text(.8,i-.2,sprintf('%s\n',textbank{5+i}{:}),'HorizontalAlignment','right','FontSize',14);
end
legend(ActivityLabels(1:5),'Location','Best');

SUPERPLOT(1,4,1,3,.2,.1,.05);  set(gca,'FontSize',14)
for i = 1:10
    h(i) = plot(nanmean(participantconf(:,i,:),3),'LineWidth',1,'Color',cmp(i,:)); hold on;
end
ylim([-4.1 -.7]); xlim([1 size(participantconf,1)]); title('Participant');
plot(nanmean(nanmean(participantconf,3),2),'LineWidth',3,'Color','k');
set(gca,'Yticklabel',[],'FontSize',12)
xlabel('Session Number','FontSize',16);
legend(h(6:10),ActivityLabels(6:10),'Location','Best');

SUPERPLOT(1,12,1,10,.2,.1,.05); 
scatter(ones(10,1)+.02*rand(10,1),[nanmean(nanmean(researcherconf,3))]',[],cmp,'filled'); hold on;
scatter(1.1*ones(10,1)+.02*rand(10,1), [nanmean(nanmean(participantconf,3))]',[],cmp,'filled');
set(gca,'Yticklabel',[],'Xticklabel',[],'FontSize',12)
title('Avg.')
ylim([-4.1 -.7]); xlim([.95 1.15]);
text(1,-4.3,'Res.','FontSize',12,'HorizontalAlignment','center');
text(1.1,-4.3,'Part.','FontSize',12,'HorizontalAlignment','center');

tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/confplot.pdf'],'pdf');

%% Log time data
tTasks = {'Turn on the Computer','Launch Adobe Connect','Calibrate Eye Tracker',...
    'Apply Cap','Apply Electrodes','OpenBCI2000','Check Impedances',...
    'Perform a P300 Run','Update Classifier','Clean System'}
clear plogtimes rlogtimes
for tt = 1:length(P.code)
    %Load log from each session
    allsess = dir(['data/' P.code{tt} '*']);
    allsessRCAP = [P.initial{tt}; P.tbcisess{tt}; P.interaction{tt}(~isnan(P.interaction{tt}))];
    plogtimes{tt} = [];
    for sess = 1:length(allsess)
        tmpsess = allsess(sess).name;
        fnme = dir(['data/' tmpsess '/LogData*']);
        fpp = fopen(['data/' tmpsess '/' fnme.name],'r');
        if fpp~=-1
            formatSpec = '%s%[^\n\r]';
            dataArray = textscan(fpp, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
            txtt{tt,sess} = dataArray{1};
            if ~isempty(find(strcmp(txtt{tt,sess},'SYSTEM RESET')))
                disp(['SYSTEM RESET found in ' tmpsess]);
            end
                
            fclose(fpp);
            for i = 1:length(tTasks)
                tmpi = find(strcmp(txtt{tt,sess},tTasks{i}));
                
                if isempty(tmpi)
                    plogtimes{tt}(sess,i) = NaN;
                else
                    plogtimes{tt}(sess,i) = sum(cellfun(@(x) str2num(x), txtt{tt,sess}(tmpi+1)))/60;
                end
                
            end
        else
            disp(['NO DATA EXISTS, ' P.code{tt} 'S' num2str(sess)]);
            plogtimes{tt}(sess,:) = NaN(1,10);
        end
        
        
    end
    rlogtimes{tt}= [cell2mat(PData(allsessRCAP,IND.rlog_time));
        NaN(10-length(allsessRCAP),length(IND.rlog_time))];
    plogtimes{tt} = [plogtimes{tt};
        NaN(10-sess,length(IND.rlog_time))];
end

rlogtimes_ALL = cat(3,rlogtimes{:});
plogtimes_ALL = cat(3,plogtimes{:});

hh=figure('Position',[100 100 800 500]);
barweb_AG(nanmean(rlogtimes_ALL,3),nanstd(rlogtimes_ALL,[],3),[],...
    {'I','S1','S2','S3','S4','S5','S6','S7','S8','Int'});
ylim([0 20]); xlim([.5 size(rlogtimes_ALL,1)+.5]); set(gca,'FontSize',14);
xlabel('Session'); ylabel('Time (min.)'); title('Researcher Log');
legend(ActivityLabels,'Location','EastOutside')
tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/rtimelog.pdf'],'pdf');


hh=figure('Position',[100 100 800 500]);
barweb_AG(nanmean(plogtimes_ALL,3),nanstd(plogtimes_ALL,[],3),[],...
    {'I','S1','S2','S3','S4','S5','S6','S7','S8','Int'}); 
ylim([0 20]); xlim([.5 size(plogtimes_ALL,1)+.5]); set(gca,'FontSize',14);
xlabel('Session'); ylabel('Time (min.)'); title('Guide Log');
legend(ActivityLabels,'Location','EastOutside')
tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/ptimelog.pdf'],'pdf');

%% Time to first run analysis
%This is steps 1:7 in the logs
%Only calculate for telebci sessions (2:9)

rFRA = squeeze(nansum(rlogtimes_ALL(2:9,1:7,:),2));
rFRA(rFRA==0) = NaN;
pFRA = squeeze(nansum(plogtimes_ALL(2:9,1:7,:),2));
pFRA(pFRA==0) = NaN;
mxx = max(max([rFRA;pFRA]));

xlab = 1:8;
xlabm = repmat(xlab',1,size(pFRA,2));
xlabr = xlabm+.15-.3*rand(size(pFRA));
cmp = parula(size(pFRA,2));
plab = repmat(1:size(pFRA,2),8,1);
hh = figure;
subplot(121); hold on;
remind = ~isnan(pFRA(:));
scatter(xlabr(remind),pFRA(remind),300,cmp(plab(remind),:),'Marker','.');
title(sprintf('%s\n','Log-reported','Time to first P300 run'));
ylim([0 1.1*mxx]); xlim([.5 8.5]);
%%TRENDLINE (Quadratic?)
p = polyfit(xlabm(remind),pFRA(remind),2)
x1 = linspace(1,8,100);
y1 = polyval(p,x1);
plot(x1,y1)
subplot(122);  hold on;
remind = ~isnan(rFRA(:));
scatter(xlabr(remind),rFRA(remind),300,cmp(plab(remind),:),'Marker','.');
title(sprintf('%s\n','Researcher-reported','Time to first P300 run'));
ylim([0 1.1*mxx]); xlim([.5 8.5]);
%%TRENDLINE (Quadratic?)
p = polyfit(xlabm(remind),rFRA(remind),2)
x1 = linspace(1,8,100);
y1 = polyval(p,x1);
plot(x1,y1)


% % Tried fitting to an exponential as well.  Did not understand why to use
% % this function -- gets similar results as quadratic
% % https://www.mathworks.com/help/optim/examples/nonlinear-data-fitting.html
% F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);
% x0 = [1 1 1 0];
% [x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,xlabm(remind),rFRA(remind))
% plot(xlab,F(x,xlab))






tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/time2firstrun.pdf'],'pdf');


%% Impedance data
ElecNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};

clear IMP
for tt = 1:length(P.code)
    %Load log from each session
    allsess = dir(['data/' P.code{tt} '*']);

    for sess = 1:length(allsess)
        tmpsess = allsess(sess).name
        fnme = dir(['data/' tmpsess '/impedance*']);
        if isempty(fnme)
            fpp=-1;
        else
            fpp = fopen(['data/' tmpsess '/' fnme(end).name],'r')
        end
        if fpp~=-1
            formatSpec = '%s%s%s%s%s%s%s%s%[^\n\r]';
            dataArray = textscan(fpp, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
            imptmp = dataArray{1};
            fclose(fpp)
            
            
            clear e_ind e_loc
            for jj = 1:length(ElecNames)
                for ll = 1:length(  imptmp)
                    e_ind(jj,ll) = ~isempty(strfind(imptmp{ll},ElecNames{jj}));
                end
                e_loc{jj} = find(e_ind(jj,:)==1);
            end
            e_loc_end = cellfun(@(x) x(end),e_loc);
            IMPS = imptmp(e_loc_end);
            expr = [':' '(.*?)' 'kOhm'];
            tok = cellfun(@(x) regexp(x,expr,'tokens'),IMPS);
            IMP(tt,sess,:)=cell2mat(cellfun(@(x) str2num(cell2mat(x)),tok,'UniformOutput',false));
          
            
            %         [CData, NData, Spec] =  ShowMeTheP300(tempName,tempSession{sess},...
            %             tempRun(ii),0,2,0, full_code);
            %
            %         CD(kk,:,:) = mean(CData,3);
            %         ND(kk,:,:) = mean(NData,3);
            %         SD(kk,:,:) = Spec;
           
        else
            IMP(tt,sess,:) = NaN(1,length(ElecNames));
        end
        
    end
    labelP{tt} = [P.code{tt} sess];
end

hh=figure('Position',[100 100 1200 400]);
for i = 1:length(P.code)
    SUPERPLOT(1,length(P.code),1,i,.17 ,.15 ,.15)
    imagesc(squeeze(IMP(i,:,:))',[0 50]);
    set(gca,'Yticklabel',[],'Fontsize',16);
    if i == 1
        xlabel('Session Number');
        for j = 1:8
            text(0.4,j,ElecNames{j},'HorizontalAlignment','right','Fontsize',16);
        end
    end
end
cc = colorbar; ylabel(cc, 'Impedance (k\Omega)')
tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/impedances.pdf'],'pdf');

%% Gaze Data

gazetracking(P.code)

%% P300 data
clear P300Fisher P300avg_s Navg_s P3_GA P3_GS...
    N_GA N_GS P3_GA_s N_GA_s P3_GA_s2 N_GA_s2
close all

fltext = '';
fs = 250;
% drange = (round(.45*fs):round(1.2*fs)); %This is the window that is in
% the data
drange = fix(1:ceil(.75*fs))/fs;
ChanNames = {'Fz','Cz','P3','Pz','P4','PO7','PO8','Oz'};



for sj = 1
    if sj == 1
        ppl = P.code; pmark = '';
    elseif sj == 2
        ppl = ppl_S2(crng_S2); pmark = 'c';
    end
   
    aa = 1;
    cfm = [];
    for REF = [0 1]
        clear C_K N_K P300avg Navg P300Quality
        for ii = 1:length(ppl)
            try
                load(['results\'...
                    ppl{ii} fltext '_P300feats_Ref' num2str(REF)]);
                
                C_K(aa,ii) = C_keep;
                N_K(aa,ii) = N_keep;
                   %         end
            catch
                disp('eh')
                cfm = NaN(length(drange),length( ChanNames));
                nfm = NaN(length(drange),length( ChanNames));
                cfs = NaN(length(drange),length( ChanNames));
                nfs = NaN(length(drange),length( ChanNames));
                C_K(aa,ii) = NaN;
                N_K(aa,ii) = NaN;
            end
            P300avg(aa,ii,:,:) = cfm';
            Navg(aa,ii,:,:) = nfm';
            P300Quality(aa,ii,:,:) = (cfm'-nfm').^2./(cfs'+nfs');
            
            %Shouldnt try to do Fisher scoring unless i have individual
            %trials
            %             FSdata = cat(1,C_m(:,:,jj),N_m(:,:,jj));
            %             FSlabels = [ones(size(C_m,1),1); 2*ones(size(N_m,1),1)];
            %             temp = fsFisher(FSdata,FSlabels);
            %             P300Fisher(aa,ii,jj,:) = temp.W;
            
            
            
        end
        
        
        %     %Normalize P300 for each subject between 0 and 1
        %     temp = squeeze(P300avg(aa,:,:,:));
        %     temp2 = temp-repmat(min(min(temp,[],3),[],1),[size(temp,1),1,256]);
        %     P300avg_s(aa,:,:,:) = temp2./repmat(max(max(temp2,[],3),[],1),[size(temp2,1),1,256]);
        %     temp3 = squeeze(Navg(aa,:,:,:));
        %     temp4 = temp3-repmat(min(min(temp,[],3),[],1),[size(temp,1),1,256]);
        %     Navg_s(aa,:,:,:) = temp4./repmat(max(max(temp2,[],3),[],1),[size(temp2,1),1,256]);
        %
        %     %Normalized Grand Averages
        %     P3_GA_s2(aa,:,:) = squeeze(mean(P300avg_s(aa,pplrng,:,:),2));
        %     N_GA_s2(aa,:,:) = squeeze(mean(Navg_s(aa,pplrng,:,:),2));
        
        
        aa = aa+1;

    
    if size(cfm,2)==22
        psz1 = 5; psz2=5;
        eleclocs = [2 4 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 24];
        elecplace = [mod(eleclocs,5); floor((100-eleclocs)/5)];
        elecplace(elecplace==0)=5;
        classic_chan = [10 15];
        ChP = [21 5 10 15];
    elseif size(cfm,2)==9
        psz1 = 3; psz2 = 2;
        eleclocs = 1:6;
        elecplace = [1 2 1 2 1 2;3 3 2 2 1 1];
        ChanNames = {'F3','F4','C3','C4','P3','P4','lEOG','cEOG','rEOG'};
        classic_chan = 1:6;
    elseif size(cfm,2)==5
        psz1 = 1; psz2 = 2;
        eleclocs = 1:2;
        elecplace = [1 2;1 1];
        ChanNames = {'CSP1','CSP2','lEOG','cEOG','rEOG'};
    elseif size(cfm,2)==8
        psz1 = 3; psz2 = 4;
        eleclocs = [2 5 7 8 9 10 11 12];
        elecplace = [2 2 1 2 3 1 2 3;1 2 3 3 3 4 4 4]; 
        ChP = [1 2 4 8]
    end
    
    
    %Plot P300 averages for each person for
    %channel 21 cEOG and channel 10 (Cz) and channel 18 (O1)
    
    %Only need to do this when there is only one participant
%     P300avg = repmat(P300avg,1,3,1,1);
%    Navg = repmat(Navg,1,3,1,1);
%    ppl = repmat(ppl,3,1);
%    C_K = repmat(C_K,1,3);
%    N_K = repmat(N_K,1,3);
   
   
   cmp = parula(100);
    ALS_plottrace([ pmark 'P300_Art0_Ref' num2str(REF)],...
        {squeeze(P300avg(1,:,:,:)), squeeze(Navg(1,:,:,:))},...
        {C_K(1,:), N_K(1,:)},1,ppl,1:length(ppl), {[0 0 0],...
        [.6 .6 .6]},1,[10 10 10 10],...
        ChP,ChanNames(ChP),drange,max([length(P.code)])+1,0,[1 2]);
    
        ALS_plottrace([ pmark 'Qual_Art0_Ref' num2str(REF)],...
        {squeeze(P300Quality(1,:,:,:))},...
        {C_K(1,:)},1,ppl,1:length(ppl), {[0 0 0]},0,[1 1 1 1],...
        ChP,ChanNames(ChP),drange,max([length(P.code)])+1,0,[1 2]);
    
    aa = 1;
    end
end

%% Communication type over sessions

ctmp = cell2mat(PData(:,IND.log_comm));
ctmpr = cell2mat(PData(:,IND.rlog_comm));
fullcommr = []; fullcomm = [];
for i = 1:length(P.code)
    tmp = find(strcmp(P.fullcode,P.code{i}) & ...
        cellfun(@(x) ~isempty(x),regexp(PData(:,2),'init'))+cellfun(@(x) ~isempty(x),regexp(PData(:,2),'tele')));
tmp2 = ctmp(tmp,:);
fullcomm = cat(3,fullcomm,[tmp2; NaN(10-size(tmp2,1),size(tmp2,2))])
tmp3 = ctmpr(tmp,:);
fullcommr = cat(3,fullcommr,[tmp3; NaN(10-size(tmp3,1),size(tmp3,2))])
end

xlab = repmat(1:10,3,1)+repmat([-.2; 0; .2],1,10)
cmp = parula(3);
hh = figure('Name','Reported interaction type per session','Position',[100 100 700 200])
subplot(121); hold on;
for i = 1:3
eval(['h' num2str(i)]) = bar(xlab(i,:),(nansum(fullcommr(:,i,:),3)./sum(~isnan(fullcommr(:,i,:)),3)),.4/2);
set(eval(['h' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend({'Online','In-person','Phone'}); xlim([0.5 9.5]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','T1','T2','T3','T4','T5','T6','T7','T8'},...
    'HorizontalAlignment','center')
title('Researcher'); ylabel(sprintf('%s\n','Fraction of participants','having interaction'));
subplot(122); hold on;
for i = 1:3
eval(['h' num2str(i)]) = bar(xlab(i,:),(nansum(fullcomm(:,i,:),3)./sum(~isnan(fullcomm(:,i,:)),3)),.4/2);
set(eval(['h' num2str(i)]),'FaceColor',cmp(i,:)); 
end
title('Participant Team');  xlim([0.5 9.5]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','T1','T2','T3','T4','T5','T6','T7','T8'},...
    'HorizontalAlignment','center')

tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/IntxType.pdf'],'pdf');

%% Help type over sessions

htmp = cell2mat(PData(:,IND.log_help(1:end-1)));
htmpr = cell2mat(PData(:,IND.rlog_help(1:end-1)));
fullhelpr = []; fullhelp = [];
for i = 1:length(P.code)
    tmp2 = htmp([P.initial{i}; P.tbcisess{i}],:);
    fullhelp = cat(3,fullhelp,[tmp2; NaN(9-size(tmp2,1),size(tmp2,2))]);
    tmp3 = htmpr([P.initial{i}; P.tbcisess{i}],:);
    fullhelpr = cat(3,fullhelpr,[tmp3; NaN(9-size(tmp3,1),size(tmp3,2))]);
end

xlab = repmat(1:9,11,1)+repmat(linspace(-.4,.4,11)',1,9);
cmp = parula(11);
gg = figure('Name','Researcher reported help type per session','Position',[100 100 700 200])
hold on;
for i = 1:11
eval(['g' num2str(i)]) = bar(xlab(i,:),(nansum(fullhelpr(:,i,:),3)./sum(~isnan(fullhelpr(:,i,:)),3)),.8/10);
set(eval(['g' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend([ActivityLabels 'Other'],'Location','eastoutside'); 
xlim([0.5 9.5]); ylim([0 1]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','S1','S2','S3','S4','S5','S6','S7','S8'},...
    'HorizontalAlignment','center')
ylabel(sprintf('%s\n','Fraction of participants','receiving help'));

tmp = get(gg,'Position');
set(gg,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(gg,['figures/ResHelpType.pdf'],'pdf');

ff = figure('Name','Participant reported help type per session','Position',[100 100 700 200])
hold on;
for i = 1:11
eval(['g' num2str(i)]) = bar(xlab(i,:),(nansum(fullhelp(:,i,:),3)./sum(~isnan(fullhelp(:,i,:)),3)),.8/10);
set(eval(['g' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend([ActivityLabels 'Other'],'Location','eastoutside'); 
xlim([0.5 9.5]); ylim([0 1]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','S1','S2','S3','S4','S5','S6','S7','S8'},...
    'HorizontalAlignment','center')
ylabel(sprintf('%s\n','Fraction of participants','receiving help'));

tmp = get(ff,'Position');
set(ff,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(ff,['figures/PatHelpType.pdf'],'pdf');



hh = figure('Name','Differences in reported help per session','Position',[100 100 700 200])
hold on;
for i = 1:11
eval(['g' num2str(i)]) = bar(xlab(i,:),((nansum(fullhelp(:,i,:),3)-nansum(fullhelpr(:,i,:),3))./sum(~isnan(fullhelp(:,i,:)),3)),.8/10);
set(eval(['g' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend([ActivityLabels 'Other'],'Location','eastoutside'); 
xlim([0.5 9.5]); ylim([-1 1]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-1.06*ones(1,9),{'I','S1','S2','S3','S4','S5','S6','S7','S8'},...
    'HorizontalAlignment','center')
ylabel(sprintf('%s\n','Difference in perception','of help recieved','(+ is dyad reporting','receiving more help)'));

tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/HelpDiff.pdf'],'pdf');



%% BCI Results

tbci = tBCIFileNames;
OnlineResults(tbci);
load('OnlineResults.mat');


tmp = cellfun(@(x) x==1, P300_r.type,'UniformOutput',false);
%Visual spelling accuracy for letters
tmp2 = cellfun(@(x,y,z) (x(z)==y(z)),P300_r.S_TEXT,P300_r.G_TEXT,tmp,'UniformOutput',false);
tmp3 = cellfun(@(x) nanmean(x),tmp2,'UniformOutput',false);
% newacc = cell2mat(tmp3);

%symbol visual speller accuracy
smp = regexp(P300_r.G_TEXT,'[A-Z][a-z]');
for n = 1:size(smp,1)
    for i= 1:size(smp,2)
        tmp4{n,i} = [];
        for j = 1:length(smp{n,i})
            tmp4{n,i}(j) = strcmp(P300_r.S_TEXT{n,i}(smp{n,i}(j):smp{n,i}(j)+2),P300_r.G_TEXT{n,i}(smp{n,i}(j):smp{n,i}(j)+2))
        end
    end
end
% newacc_s = cellfun(@(x) nanmean(x),tmp4)

%audio speller accuracy
tmp5 = cellfun(@(x) x==2, P300_r.type,'UniformOutput',false);
%Visual spelling accuracy for letters
tmp6 = cellfun(@(x,y,z) (x(z)==y(z)),P300_r.S_TEXT,P300_r.G_TEXT,tmp5,'UniformOutput',false);


tmp7 = cellfun(@(x,y,z) cat(2,x,y,z),tmp2,tmp4,tmp6,'UniformOutput',false);
newacc = cellfun(@(x) nanmean(x),tmp7);






hh = figure('Position',[100 100 300 300]) 
xx=barweb_AG(100*nanmean(newacc)',100*nanstd(newacc)'); set(gca,'FontSize',18);
xlim([.5 1.5]); ylim([0 120]);
set(gca,'Xticklabel',[]);
for i = 1:10
    text(i/12.5+.56,-2,num2str(i),'HorizontalAlignment','center')
end
xlabel('Session')
ylabel('Accuracy (%)');
set(gca,'Xticklabel',{'I','S1','S2','S3','S4','S5','S6','S7','S8','Int'})
tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/accuracy.pdf'],'pdf');


%%Alternative plot
hh = figure('Name','Scatterplot of Accuracies','Position',[100 100 600 200]);
cmp = parula(size(newacc,1));
for i = 1:size(newacc,1)
sessd = (1:10)+.2-.4*rand(1,10);
dd(i) = scatter(sessd,100*newacc(i,:),300,cmp(i,:),'Marker','.'); hold on;
line(repmat(1:10,2,1)+repmat([.2;-.2],1,10),repmat(100*nanmedian(newacc),2,1),'Color','k');
end
xlim([0.5 10.5]);
set(gca,'Xticklabel',{'I','S1','S2','S3','S4','S5','S6','S7','S8','Int'})
ylabel('Accuracy (%)');
legend(dd,P.code)
%Trendline



tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/scatter_accuracy.pdf'],'pdf');


%% What systems were used?
%Check which speller types were utilized each session.  If notepad, what 
%features were used?
ttmp = cell2mat(PData(:,[IND.log_speller IND.log_access]));
ttmpr = cell2mat(PData(:,[IND.rlog_speller IND.rlog_access]));
fulltype = []; fulltyper = [];
for i = 1:length(P.code)
    tmp = find(strcmp(P.fullcode,P.code{i}) & ...
        cellfun(@(x) ~isempty(x),regexp(PData(:,2),'init'))+...
        cellfun(@(x) ~isempty(x),regexp(PData(:,2),'tele')));
tmp2 = ttmp(tmp,:);
fulltype = cat(3,fulltype,[tmp2; NaN(9-size(tmp2,1),size(tmp2,2))]);
tmp3 = ttmpr(tmp,:);
fulltyper = cat(3,fulltyper,[tmp3; NaN(9-size(tmp3,1),size(tmp3,2))]);
end


xlab = repmat(1:9,4,1)+repmat(linspace(-.3,.3,4)',1,9)
cmp = parula(4);
gg = figure('Name','Researcher reported speller type per session','Position',[100 100 700 200])
hold on;
for i = 1:4
eval(['g' num2str(i)]) = bar(xlab(i,:),(nansum(fulltyper(:,i,:),3)./sum(~isnan(fulltyper(:,i,:)),3)),.6/3);
set(eval(['g' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend({'Trainer','Symbolic','Dual','Notepad'},'Location','eastoutside'); 
xlim([0.5 9.5]);  ylim([0 1]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','S1','S2','S3','S4','S5','S6','S7','S8'},...
    'HorizontalAlignment','center')
ylabel(sprintf('%s\n','Fraction of participants','using speller'));
tmp = get(gg,'Position');
set(gg,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(gg,['figures/ResSpType.pdf'],'pdf');

gg = figure('Name','Participant reported speller type per session','Position',[100 100 700 200])
hold on;
for i = 1:4
eval(['g' num2str(i)]) = bar(xlab(i,:),(nansum(fulltype(:,i,:),3)./sum(~isnan(fulltype(:,i,:)),3)),.6/3);
set(eval(['g' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend({'Trainer','Symbolic','Dual','Notepad'},'Location','eastoutside'); 
xlim([0.5 9.5]);  ylim([0 1]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','S1','S2','S3','S4','S5','S6','S7','S8'},...
    'HorizontalAlignment','center')
ylabel(sprintf('%s\n','Fraction of participants','using speller'));
tmp = get(gg,'Position');
set(gg,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(gg,['figures/PatSpType.pdf'],'pdf');


xlab = repmat(1:9,3,1)+repmat(linspace(-.3,.3,3)',1,9)
cmp = parula(3);
gg = figure('Name','Researcher reported notepad accessories per session','Position',[100 100 700 200])
hold on;
for i = 1:3
eval(['g' num2str(i)]) = bar(xlab(i,:),(nansum(fulltyper(:,i+4,:),3)./sum(~isnan(fulltyper(:,i+4,:)),3)),.6/2);
set(eval(['g' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend({'Text-to-Speech','Text Prediction','Prepared Text'},'Location','eastoutside'); 
xlim([0.5 9.5]);  ylim([0 1]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','S1','S2','S3','S4','S5','S6','S7','S8'},...
    'HorizontalAlignment','center')
ylabel(sprintf('%s\n','Fraction of participants','using speller'));
tmp = get(gg,'Position');
set(gg,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(gg,['figures/ResAccessores.pdf'],'pdf');

gg = figure('Name','Participant reported notepad accessories per session','Position',[100 100 700 200])
hold on;
for i = 1:3
eval(['g' num2str(i)]) = bar(xlab(i,:),(nansum(fulltype(:,i+4,:),3)./sum(~isnan(fulltype(:,i+4,:)),3)),.6/2);
set(eval(['g' num2str(i)]),'FaceColor',cmp(i,:)); 
end
legend({'Text-to-Speech','Text Prediction','Prepared Text'},'Location','eastoutside');
xlim([0.5 9.5]); ylim([0 1]);
set(gca,'Xticklabel',[]);
set(gca,'XTick',(1:8)+.5);
text([1 2 3 4 5 6 7 8 9],-.03*ones(1,9),{'I','S1','S2','S3','S4','S5','S6','S7','S8'},...
    'HorizontalAlignment','center')
ylabel(sprintf('%s\n','Fraction of participants','using speller'));
tmp = get(gg,'Position');
set(gg,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(gg,['figures/PatAccessores.pdf'],'pdf');
%% Interaction Results

intemp = cell2mat(P.interaction);
itloc = find(~isnan(intemp));
P.Ncode = PData(intemp(itloc),IND.intx_id);
P.intx_date = PData(intemp(itloc),IND.intx_date);
P.intx_conf = 5-cell2mat(PData(intemp(itloc),IND.intx_assess(1:11)));
P.intx_comf = 5-cell2mat(PData(intemp(itloc),IND.intx_assess(12:15)));
P.intx_rate = 6-cell2mat(PData(intemp(itloc),IND.intx_assess(16:21)));
P.intx_comm = 6-cell2mat(PData(intemp(itloc),IND.intx_assess(22:24)));
P.intx_answ = 6-cell2mat(PData(intemp(itloc),IND.intx_assess(25:27)));
P.intx_time = cell2mat(PData(intemp(itloc),IND.intx_assess(28:37)));
P.intx_tapp = cell2mat(PData(intemp(itloc),IND.intx_assess(38:47)));
P.intx_change = PData(intemp(itloc),IND.intx_assess(48:57));
P.intx_acc = cell2mat(PData(intemp(itloc),IND.intx_assess(58:61)));
P.intx_help = 5-cell2mat(PData(intemp(itloc),IND.intx_assess(62:65)));
P.intx_cmts = PData(intemp(itloc),IND.intx_assess(66));

P.intx_N_acc = cell2mat(PData(intemp(itloc),IND.intx_nassess(1:4)));
P.intx_N_help = 5-cell2mat(PData(intemp(itloc),IND.intx_nassess(5:8)));
P.intx_N_comm = 5-cell2mat(PData(intemp(itloc),IND.intx_nassess(9:11)));
P.intx_N_answ = 5-cell2mat(PData(intemp(itloc),IND.intx_nassess(12:14)));
P.intx_N_rate = 6-cell2mat(PData(intemp(itloc),IND.intx_nassess(15:17)));
P.intx_N_cmts = PData(intemp(itloc),IND.intx_nassess(18));

alldata = {P.intx_conf,P.intx_comf,P.intx_rate,P.intx_comm,P.intx_answ,P.intx_tapp,...
    P.intx_help,P.intx_N_comm,P.intx_N_answ,P.intx_N_rate,P.intx_N_help};
allnames = {'Participant Confidence','Participant Comfort','Participant Ratings',...
    'Patient Communicate Concerns','Patient Answer Questions','Patient Time Appropriateness',...
    'Patient Helpful Accessories',...
    'Nurse Communicate Concerns','Nurse Answer Questions','Nurse Ratings','Nurse Helpful Accessories'};
alllabels = {{'VC','SC','LC','NC'},{'VC','C','U','VU'},{'VH','H','M','L','VL'},...
    {'VW','W','P','VP','DK'},{'VW','W','P','VP','DK'},{'WL','L','A','S','WS'},...
    {'VH','H','U','VU'},...
    {'VW','W','P','VP','DK'},{'VW','W','P','VP','DK'},{'VH','H','M','L','VL'},{'VH','H','U','VU'}};
allitems = {{'Comp','Adobe Con','Eye','BCI2000','Cap','Elec','Imp','P300','Class','Remove','Care'},...
    {'Cap','Elec','Eye','P300'},{'Efficacy','Accuracy','Speed','Setup','Use','Portability'},...
    {'w/o AAC','w/AAC','w/P300'},{'w/o AAC','w/AAC','w/P300'},{'Comp','Adobe Con','Eye','BCI2000','Cap','Elec','Imp','P300','Class','Remove','Care'},...
    {'Text-to-Speech','Word Prediction','Symbolic Communication','Prepared Text'},...
    {'w/o AAC','w/AAC','w/P300'},{'w/o AAC','w/AAC','w/P300'},...
    {'Efficacy','Accuracy','Speed'},{'Text-to-Speech','Word Prediction','Symbolic Communication','Prepared Text'}};
allrange = {[4 3 2 1],[4 3 2 1],[5 4 3 2 1],[5 4 3 2 1],[5 4 3 2 1],[5 4 3 2 1],...
    [4 3 2 1],[5 4 3 2 1],[5 4 3 2 1],[5 4 3 2 1],[4 3 2 1]};

for i = 1:length(alldata)
    hh = figure('Name',allnames{i},'Position',[100 100 300 200])
    tmpd = alldata{i};
    tmpn = allnames{i};
    tmpl = alllabels{i};
    tmpx = repmat(1:size(tmpd,2),size(tmpd,1),1)-.2*rand(size(tmpd))
    cmp = parula(size(tmpd,2));
    cmpfull = [];
    for j = 1:size(cmp,1)
        cmpfull = [cmpfull; repmat(cmp(j,:),size(tmpd,1),1)];
    end
    scatter(tmpx(:),tmpd(:),100,cmpfull,'.');
    xx = line(repmat(1:size(tmpd,2),2,1)+repmat([-.2;.2],1,size(tmpd,2)),...
        repmat(nanmedian(tmpd),2,1))
    for j = 1:length(xx)
        xx(j).Color = cmp(j,:)
    end
    set(gca,'Xtick',(1:size(tmpd,2)-1)+.5,'Xticklabels',[],'Yticklabels',[]);
    
    text(1:length(allitems{i}),min(allrange{i}-.6)*ones(length(allitems{i}),1),...
        allitems{i},'Rotation',-45); %xlabels
    ylim([min(allrange{i})-.5 max(allrange{i})+.5]);
    xlm = [.5 size(tmpd,2)+.5]
    xlim(xlm);
    diff(xlm)
    text(.5*ones(length(allrange{i}),1)-diff(xlm)/30,allrange{i},alllabels{i},...
        'HorizontalAlignment','right'); %ylabels
    
    tmp = get(hh,'Position');
    set(hh,...
        'DefaultAxesFontSize',5,...
        'DefaultLineLineWidth',2,...
        'PaperUnits','points',...
        'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
        'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
    saveas(hh,['figures/' allnames{i} '.pdf'],'pdf');
    
end

    



