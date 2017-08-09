%Online Results - average all  runs (allruns = 1), or do selective
%averaging (online last two runs of session, allruns = 0), or just take the
%best run (allruns = 2, NOT COMPLETE YET)


%Can just run this file after running
%[P300 MI] = ALSStudyFileNames()

%This only works currently for allruns = 1
function [P300_r] = OnlineResults(P300)

runloc = cd;
[~, sl] = regexp(runloc,'\Users\');
se = regexp(runloc,'Box Sync');
usrnm = runloc(sl+1:se-2);
runloc = runloc(1);


%Check for current file
try
    load(['OnlineResults_ ' num2str(allruns) '.mat']);
    if P300_r.fls == length(P300.Names)
        disp('Online results up to date -- not running');
        return
    end
catch
end

%% P300 Online results

%Easiest way is to access the logfiles

nms = P300.Names;
P300_r.Subjs = nms;
P300_r.fls = length(P300.Names);
ALLTEXT = []; ALLTYPE = []; ALLSUBJ = [];
P300_r.sacc = cell(length(nms),10);
P300_r.G_TEXT = cell(length(nms),10);
P300_r.S_TEXT = cell(length(nms),10); 
P300_r.type = cell(length(nms),10);

for ss = 1:length(nms)
    nms{ss}
    
    
    %s_indices =  P300.Feedback{ss} == 1;
    
    for ssn = 1:length(P300.Session{ss})
        ssn
        allind = P300.Run{ss}{ssn}
        
        for rr = 1:length(allind)
            
            fid = fopen(['data\' nms{ss} P300.Session{ss}{ssn} '\'...
                nms{ss} 'S' P300.Session{ss}{ssn} 'R' num2str(str2num(allind{rr})) '+_mylogfile.txt'],'r');
            if fid==-1
                disp('no log file -- checking if audio speller');
                
                [a,b,c,d] = load_bcidat(['data\' nms{ss} P300.Session{ss}{ssn} '\'...
                    nms{ss} 'S' P300.Session{ss}{ssn} 'R' allind{rr}],'-calibrated');
                %remove 8 seconds
                remind = 1:8*c.SamplingRate.NumericValue;
                b.StimulusCode(remind) = 0;
                b.StimulusType(remind) = 0;
                b.SelectedStimulus(remind) = 0;
                if size(c.Stimuli.Value,2)==8
                    disp('audio speller');
                else
                    disp('ERROR -- visual nor audio speller recognized')
                    return;
                end
                StimCode = b.StimulusCode(b.StimulusType>0);
                G_t = StimCode(1:c.NumberOfSequences.NumericValue*...
                    c.StimulusDuration.NumericValue/1000*c.SamplingRate.NumericValue:end)';
                ResCode = b.SelectedStimulus(b.SelectedStimulus>0);
                S_t = ResCode(1:875:end)'; % Not sure where the 4s - 500ms comes from
                
                audiocodes = {'y','n','e','a','c','f','b','d'};
                if isempty(G_t)
                    disp('Free Spelling run... abort')
                else
                    P300_r.sacc{ss,ssn} = [P300_r.sacc{ss,ssn} G_t==S_t];
                    P300_r.type{ss,ssn} = [P300_r.type{ss,ssn} 2*ones(1,length(S_t))];
                    disp(['Goal text = ' cat(2,c.Stimuli.Value{1, G_t})...
                        ', Selected Text = ' cat(2,c.Stimuli.Value{1, S_t})])
                    P300_r.G_TEXT{ss,ssn} = [P300_r.G_TEXT{ss,ssn} cat(2,audiocodes{G_t})];
                    P300_r.S_TEXT{ss,ssn} = [P300_r.S_TEXT{ss,ssn} cat(2,audiocodes{S_t})];
                    ALLTEXT = [ALLTEXT [cat(2,audiocodes{G_t});cat(2,audiocodes{S_t})]];
                    ALLSUBJ = [ALLSUBJ repmat(ss,1,length(S_t))];
                end
                
                   
                   
                   
            else  %Proceed with calculating online results from checkerboard speller
                itt=textscan(fid,'%s','delimiter','\n');
                fclose(fid)
                itt = itt{:};
                clear G_i S_i
                for text_line = 1:length(itt)
                    G_i(text_line) = ~isempty(strfind(itt{text_line},'Goal Text,'));
                    S_i(text_line) = ~isempty(strfind(itt{text_line},'Selected Text'));
                end
                if sum(G_i)==0
                    disp('Free Spelling run... abort')
                else
                    G_i = find(G_i==1);
                    G_t = regexp(itt{G_i(end)},'=');
                    G_t = itt{G_i(end)}(G_t+2:end);
                    
                    S_t = itt(S_i);
                    S_t = cellfun(@(x) x(18:end),S_t,'UniformOutput',false)';
                    if ss == 1 && ssn == 8 && rr == 2
                        S_t(4) = []; %Manually remove the backspace from this trial
                    end
                    
                    if isempty(regexp(G_t,'[a-z]')) %Letter Speller
                    S_t = cell2mat(S_t);
                    P300_r.sacc{ss,ssn} = [P300_r.sacc{ss,ssn} G_t==S_t];
                    P300_r.type{ss,ssn} = [P300_r.type{ss,ssn} ones(1,length(S_t))];
                    else %Symbolic Speller
                    G_t = strsplit(G_t);  
                    G_t = G_t(1:end-1);
                    G_t = cellfun(@(x) [x ' '], G_t,'UniformOutput',false); %Pad symbols of length two to three
                    %Just incase there is a ^ or $ or! spelled, make the
                    %symbol three times longer
                    S_t = cellfun(@(x) repmat(x,1,3),S_t,'UniformOutput',false);
                    S_t = cellfun(@(x) x(1:3),S_t,'UniformOutput',false)
                    G_t = cellfun(@(x) x(1:3),G_t,'UniformOutput',false)
                    P300_r.sacc{ss,ssn} = [P300_r.sacc{ss,ssn} cellfun(@(x,y) strcmp(x,y),G_t,S_t)];
                    P300_r.type{ss,ssn} = [P300_r.type{ss,ssn} 2*ones(1,length(S_t))]
                    G_t = cat(2,G_t{:});
                    S_t = cat(2,S_t{:});
                    end  
                    
                 
                    disp(['Goal text = ' G_t ', Selected Text = ' S_t])
                    
                   
                    P300_r.G_TEXT{ss,ssn} = [P300_r.G_TEXT{ss,ssn} G_t];
                    P300_r.S_TEXT{ss,ssn} = [P300_r.S_TEXT{ss,ssn} S_t];
                    ALLTEXT = [ALLTEXT [G_t;S_t]];
                    ALLSUBJ = [ALLSUBJ repmat(ss,1,length(S_t))];
                end
            end
            
        end
        
    end
end

        
%         s_indices =  logical(cell2mat(P300.Feedback{ss}{ssn})) & ....
%             ~logical(cell2mat(P300.FreeSpell{ss}{ssn}));
%         s_type = cell2mat(P300.Type{ss}{ssn}(s_indices));
%         s_run = P300.Run{ss}{ssn}(s_indices);
%         
%         %     if allruns == 0
%         %         s_s = unique(s_sess);
%         %         s_rind = [];
%         %         for uu = 1:length(s_s)
%         %             temp = find(strcmp(s_sess,s_s(uu))==1);
%         %             if length(temp)>2
%         %                 s_rind = [s_rind temp(end-1:end)];
%         %             else
%         %                 s_rind = [s_rind temp];
%         %             end
%         %         end
%         %     else
%         s_rind = 1:length(s_run);
%         %     end
%         
%         
%         for si = 1:length(s_rind)
%             si
%             if runloc == 'C'
%                 if strcmp(usrnm,'ageronimo')
%                 fid = fopen(['C:\Users\ageronimo\Documents\BCI2000_305_VS2012\'...
%                     'data\' nms{ss} P300.Session{ss}{ssn} '\'...
%                     nms{ss} 'S' P300.Session{ss}{ssn} 'R' num2str(str2num(s_run{si})) '+_mylogfile.txt'],'r');
%                 elseif strcmp(usrnm,'BCIcomputer01')
%                      fid = fopen(['C:\Users\BCIcomputer01\Documents\40647StudyData\'...
%                     nms{ss} P300.Session{ss}{ssn} '\'...
%                     nms{ss} 'S' P300.Session{ss}{ssn} 'R' num2str(str2num(s_run{si})) '+_mylogfile.txt'],'r');
%                 end 
%                 if fid==-1
%                     fid = fopen(['P:\ALS Proj Data\' nms{ss} P300.Session{ss}{ssn} '\'...
%                         nms{ss} 'S' P300.Session{ss}{ssn} 'R' num2str(str2num(s_run{si})) '+_mylogfile.txt'],'r');
%                 end
%             else
%                 fid = fopen(['/gpfs/work/a/amg5106/ALS_P300/Data/' nms{ss} P300.Session{ss}{ssn} '/'...
%                     nms{ss} 'S' P300.Session{ss}{ssn} 'R' num2str(str2num(s_run{si})) '+_mylogfile.txt'],'r');
%             end
%             if fid==-1
%                 disp('ERROR --- is the flash drive plugged in?');
%             end
%             itt=textscan(fid,'%s','delimiter','\n');
%             fclose(fid)
%             itt = itt{:};
%             clear G_i S_i
%             for text_line = 1:length(itt)
%                 G_i(text_line) = ~isempty(strfind(itt{text_line},'Goal Text,'));
%                 S_i(text_line) = ~isempty(strfind(itt{text_line},'Selected Text'));
%             end
%             G_i = find(G_i==1);
%             G_t = regexp(itt{G_i(end)},'=');
%             G_t = itt{G_i(end)}(G_t+2:G_t+5);
%             
%             S_t = itt(S_i);
%             S_t = cell2mat(cellfun(@(x) x(18),S_t,'UniformOutput',false))';
%             
%             %         disp(['Goal text = ' G_t ', Selected Text = ' S_t]);
%        
%             P300_r.sacc{ss,ssn} = [P300_r.sacc{ss,ssn} sum(G_t==S_t)/length(G_t)];
%             P300_r.G_TEXT{ss,ssn} = [P300_r.G_TEXT{ss,ssn} G_t];
%             P300_r.S_TEXT{ss,ssn} = [P300_r.S_TEXT{ss,ssn} S_t];
%             P300_r.TYPE{ss,ssn} = [P300_r.TYPE{ss,ssn} repmat(s_type(si),1,length(S_t))];
%             ALLTEXT = [ALLTEXT [G_t;S_t]];
%             ALLTYPE = [ALLTYPE repmat(s_type(si),1,length(S_t))];
%             ALLSUBJ = [ALLSUBJ repmat(ss,1,length(S_t))];
%             hi=2;
%         end
%         
%     end
%     P300_r.n{ss} = [length(P300_r.G_TEXT{ss}) length(P300_r.S_TEXT{ss})];
%     if diff(P300_r.n{ss})~=0
%         disp(['Difference in length of goal and selected text for subj' nms{ss}]);
%         return
%     end
%     
%     P300_r.acc{ss} = sum(P300_r.G_TEXT{ss}==P300_r.S_TEXT{ss})/...
%         length(P300_r.G_TEXT{ss});
% end
% P300_r.alltext = ALLTEXT;
% P300_r.alltype = ALLTYPE;
% P300_r.allsubj = ALLSUBJ;
% 
% %% Evaluate trends in online accuracy
% code = ['A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';...
%     'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z';'_';'0';'1';'2';'3';'4'];
% locs = [1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6 6;...
%     1 2 3 4 5 6 1 2 3 4 5 6 1 2 5 6 1 2 5 6 1 2 3 4 5 6 1 2 3 4 5 6];
% p3type = {'FS','RS','CS','FF'};
% randtype = [1/32,1/8,1/32,1/32];
% data2wayA = [];
% factor12wayA = [];
% factor22wayA = [];
% data2wayB = [];
% factor12wayB = [];
% factor22wayB = [];
% for jj = 1:length(p3type)
%     tempind = find(ALLTYPE==jj)
%     temptext = ALLTEXT(1,tempind)
%     %Were certain letters easier to spell? (P300AccByLetter)
%     for i = 1:length(code)
%         temp = strfind(ALLTEXT(1,tempind),code(i));
%         n_code(i) = length(temp);
%         acc_code(i) = sum(ALLTEXT(2,tempind(temp))==code(i))/length(temp);
%         
%         acc_matrix(locs(1,i),locs(2,i))=acc_code(i);
%         %         %Only keep codes which were used 10+ times
%         %         if n_code(i)<10
%         %             acc_matrix10(locs(1,i),locs(2,i))=0;
%         %         else
%         %             acc_matrix10(locs(1,i),locs(2,i))=acc_code(i);
%         %         end
%     end
%     
%     nn = figure('Name',[p3type{jj} 'Accuracy per letter']);
%     imagesc(acc_matrix,[0 1])
%     text(locs(2,:)',locs(1,:)',code,'FontSize',20,'HorizontalAlignment','center')
%     text(locs(2,:)',locs(1,:)+.3',num2str(n_code'),'FontSize',10,'HorizontalAlignment','center')
%     colormap(gray); tt = colorbar; ylabel(tt,'Accuracy');
%     set(gca,'Xticklabel',[],'Yticklabel',[]);
%     
%     tmp = get(nn,'Position');
%     set(nn,...
%         'DefaultAxesFontSize',5,...
%         'DefaultLineLineWidth',2,...
%         'PaperUnits','points',...
%         'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
%         'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
%     saveas(nn,['FeatureResults/figs/' p3type{jj} 'Accuracy per letter' '.pdf'],'pdf');
%     
%     %Was there an improvement over sessions? (P300AccBySession)
%     s_mean = cellfun(@(x,y,z) mean(x(z==jj)==y(z==jj)), P300_r.G_TEXT,...
%         P300_r.S_TEXT,P300_r.TYPE,'UniformOutput',false);
%     emptycels = (cellfun(@(x,y,z) isempty(mean(x(z==jj)==y(z==jj))), P300_r.G_TEXT,...
%         P300_r.S_TEXT,P300_r.TYPE,'UniformOutput',false));
%     s_mean(cell2mat(emptycels))={NaN};
%     s_mean = cell2mat(s_mean);
%     s_mean_less = s_mean(sum(s_mean,2)~=0,:);
%     s_subjs_less = P300_r.Subjs(sum(s_mean,2)~=0);
%     clrmp = jet(size(s_mean_less,1))
%     pat_ind = cellfun(@(x) strcmp(x(1),'P'),s_subjs_less);
%     s_mean_pat = nanmean(s_mean_less(pat_ind,:));
%     s_err_pat = nanstd(s_mean_less(pat_ind,:))/sqrt(sum(pat_ind));
%     con_ind = cellfun(@(x) strcmp(x(1),'C'),s_subjs_less);
%     s_mean_con = nanmean(s_mean_less(con_ind,:));
%     s_err_con = nanstd(s_mean_less(con_ind,:))/sqrt(sum(con_ind));
%     
%     %% Old -- lineplot
% %     mm = figure('Position',[100 100 400 750],'Name',[p3type{jj} 'Accuracy by subject']);
% %     hold on;
% %     set(gca, 'ColorOrder', clrmp);
% %     plot(s_mean_less'); hold on;
% %     plot(s_mean_pat,'b','LineWidth',2);
% %     plot(s_mean_con,'r','LineWidth',2);
% %     plot([1,1]*randtype(jj),'--k'); ylim([0 1]); xlim([0 3]);
% %     xlabel('P300 session');
% %     ylabel('P300 online accuracy');
% %     legend([nms(sum(s_mean,2)~=0) 'P_{Avg}','C_{Avg}'],'location','northeastoutside')
% %     
% %     tmp = get(mm,'Position');
% %     set(mm,...
% %         'DefaultAxesFontSize',5,...
% %         'DefaultLineLineWidth',2,...
% %         'PaperUnits','points',...
% %         'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
% %         'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
% %     saveas(mm,['FeatureResults/P300/' p3type{jj} 'Accuracy by subject' '.pdf'],'pdf');
%     
% 
% %DO ANOVA HERE TO LOOK FOR GROUP/SESSION FACTORS
% data = s_mean_less(:);
% factor1 = repmat([1 2],size(s_mean_less,1),1); factor1 = factor1(:); %sess1/2
% factor2 = repmat(pat_ind'+1,1,size(s_mean_less,2)); factor2 = factor2(:); %pat/cont
% [p,atab,stats,terms] = anovan(data,{factor1 factor2},'varnames',...
%     {'Session','Group'},'model','interaction','display','off')
% 
% 
% mm = figure('Position',[100 100 300 400],'Name',[p3type{jj} 'Accuracy by subject']);
% barweb_AG(100*[s_mean_pat; s_mean_con],100*[s_err_pat; s_err_con],.9,{'Patients','Controls'},[p3type{jj}],[],[],[.3 .3 .3;.7 .7 .7]);    hold on;
% plot([0 3],[100,100]*randtype(jj),'--k'); ylim([0 109]); xlim([.5 2.5]);
% if p(1)<.05
%     line([.75  1.05; .95 1.25],100*repmat(s_mean_pat+s_err_pat,2,1)*1.03.*ones(2),'Color','k');
%     line([.85  .85],100*[(s_mean_pat(1)+s_err_pat(1))*1.03 .75],'Color','k');
%     line([1.15  1.15],100*[(s_mean_pat(2)+s_err_pat(2))*1.03 .75],'Color','k');
%     line([.85  1.15],100*[.75 .75],'Color','k');
%     line([1.75  2.05; 1.95 2.25],100*repmat(s_mean_con+s_err_con,2,1)*1.03.*ones(2),'Color','k');
%     line([1.85  1.85],100*[(s_mean_con(1)+s_err_con(1))*1.03 .96],'Color','k');
%     line([2.15  2.15],100*[(s_mean_con(2)+s_err_con(2))*1.03 .96],'Color','k');
%     line([1.85  2.15],100*[.96 .96],'Color','k');
%     if p(1)<.01
%         text([1 2],[76 97],'**','HorizontalAlignment','center')
%     else
%         text([1 2],[76 97],'*','HorizontalAlignment','center')
%     end 
% end
% 
% if p(2)<.05
%       line([.75  1.75; 1.25 2.25],[80 101; 80 101],'Color','k');
%     line([1 1],[80 105],'Color','k');
%     line([2 2],[101 105],'Color','k');
%     line([1 2],[105 105],'Color','k');
%    if p(1)<.01
%         text(1.5,106,'**','HorizontalAlignment','center')
%     else
%         text(1.5,106,'*','HorizontalAlignment','center')
%    end 
% end
%     
%  
% % mm = figure('Position',[100 100 300 400],'Name',[p3type{jj} 'Accuracy by groupsession']); hold on;
% % plot(s_mean_less(pat_ind,:)','Color',[.5 .5 1],'Marker','o');
% % plot(nanmean(s_mean_less(pat_ind,:)),'b','LineWidth',2)
% % plot(s_mean_less(con_ind,:)','Color',[1 .5 .5]);
% % plot(nanmean(s_mean_less(con_ind,:)),'r','LineWidth',2)
% 
% 
%     if jj==1   
%     ylabel('Online accuracy (%)');
%     legend('Session 1','Session 2','Location','Southeast');
%     else
%         set(gca,'Yticklabels',[]);
%     end
% 
% %     for i = 1:3
% %         text(.2,120-i*8,['p_{' atab{i+1,1} '} = ' sprintf('%.3f',atab{i+1,7})]);
% %     end
% %     
%     
%     
%     tmp = get(mm,'Position');
%     set(mm,...
%         'DefaultAxesFontSize',5,...
%         'DefaultLineLineWidth',2,...
%         'PaperUnits','points',...
%         'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
%         'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
%     saveas(mm,['FeatureResults/figs/' p3type{jj} 'Accuracy by subject' '.pdf'],'pdf');
%     
%     
%   
% %%% 2 way ANOVA for Accuracy Bitrate for group and task
% 
%    s_sum = cellfun(@(x,y,z) sum(x(z==jj)==y(z==jj)), P300_r.G_TEXT,...
%         P300_r.S_TEXT,P300_r.TYPE,'UniformOutput',false);
%        s_total = cellfun(@(x,y,z) length(x(z==jj)==y(z==jj)), P300_r.G_TEXT,...
%         P300_r.S_TEXT,P300_r.TYPE,'UniformOutput',false);
%     emptycels = (cellfun(@(x,y,z) isempty(mean(x(z==jj)==y(z==jj))), P300_r.G_TEXT,...
%         P300_r.S_TEXT,P300_r.TYPE,'UniformOutput',false));
%     s_sum(cell2mat(emptycels))={NaN};
%     s_avg = sum(cell2mat(s_sum),2)./sum(cell2mat(s_total),2);
%     
%     if jj == 4
%         
%     else
%         acc_subjs_less = P300_r.Subjs(~isnan(s_avg));
%         acc_mean_less = 100*s_avg(~isnan(s_avg));
%         pat_ind = cellfun(@(x) strcmp(x(1),'P'),acc_subjs_less);
%         con_ind = cellfun(@(x) strcmp(x(1),'C'),acc_subjs_less);
%         acc_mean_pat(jj) = mean(acc_mean_less(pat_ind));
%         acc_err_pat(jj) = std(acc_mean_less(pat_ind))/sqrt(sum(pat_ind));
%         acc_mean_con(jj) = mean(acc_mean_less(con_ind));
%         acc_err_con(jj) = std(acc_mean_less(con_ind))/sqrt(sum(con_ind));
% 
%         data2wayA = [data2wayA; acc_mean_less];
%         factor12wayA = [factor12wayA; jj*ones(length(acc_mean_less),1)]; %type
%         factor22wayA = [factor22wayA; pat_ind'+1]; %pat/cont
%         
%         if jj == 2
%             bits = nansum([log2(8)+s_avg.*log2(s_avg) (1-s_avg).*log2((1-s_avg)/(8-1))],2);
%         else
%             bits = nansum([log2(32)+s_avg.*log2(s_avg) (1-s_avg).*log2((1-s_avg)/(32-1))],2);
%         end
%         bits(bits<0)=0;
%         bitrate = bits/15*60';
%         bit_mean_less = bitrate(~isnan(s_avg));
%         bit_subjs_less = P300_r.Subjs(~isnan(s_avg));
%          pat_ind = cellfun(@(x) strcmp(x(1),'P'),bit_subjs_less);
%         con_ind = cellfun(@(x) strcmp(x(1),'C'),bit_subjs_less);
%         bit_mean_pat(jj) = nanmean(bit_mean_less(pat_ind));
%         bit_err_pat(jj) = nanstd(bit_mean_less(pat_ind))/sqrt(sum(pat_ind));
%         bit_mean_con(jj) = nanmean(bit_mean_less(con_ind));
%         bit_err_con(jj) = nanstd(bit_mean_less(con_ind))/sqrt(sum(con_ind));
%         
%         data2wayB = [data2wayB; bit_mean_less];
%         factor12wayB = [factor12wayB; jj*ones(length(bit_mean_less),1)]; %type
%         factor22wayB = [factor22wayB; pat_ind'+1]; %pat/cont
%         
%     end
%     
%     
% end
% 
% 
% %Plot accuracies
% [p,atab,stats,terms] = anovan(data2wayA,{factor12wayA factor22wayA},'varnames',...
%     {'Type','Group'},'model','interaction')%,'display','off')
% [c1,m,h,nms] = multcompare(stats,'dimension',1,'display','off');
% [c2,m,h,nms] = multcompare(stats,'dimension',2,'display','off');
% 
% 
% 
% mm = figure('Position',[100 100 300 400],'Name','P300 Accuracy across groups and tasks');
% barweb_AG([acc_mean_pat; acc_mean_con],[acc_err_pat; acc_err_con],.9,{'Patients','Controls'},[],[],[],[.3 .3 .3;.6 .6 .6;.9 .9 .9]);    hold on;
% plot([0 3],max([acc_err_pat+acc_mean_pat acc_err_con+acc_mean_con]),'--k'); 
% ylim([0 1.4*max([acc_err_pat+acc_mean_pat acc_err_con+acc_mean_con])]); xlim([.5 2.5]);
% legend('FS','RS','CS','Orientation','Horizontal')
% ylabel('Accuracy (%)');
% 
% maxx = acc_mean_pat+acc_err_pat;
% maxx2 = acc_mean_con+acc_err_con;
% if p(1)<.05 %significant effect of P300 type
%     for tt=find(sign(c1(:,3).*c1(:,5))==1)'
%         idx = c1(tt,1:2)
%         line((idx-2)*.22+1,max(maxx(idx))+[.8 .8],'Color','k')
%         line(repmat((idx-2)*.22+1,2,1),max(maxx(idx))+[.5 .5;.8 .8],'Color','k')
%         text((mean(idx)-2)*.22+1,max(maxx(idx))+1.1,'*','Color','k','HorizontalAlignment','center')
%     end
%     for tt=find(sign(c1(:,3).*c1(:,5))==1)'
%         idx = c1(tt,1:2)
%         line((idx-2)*.22+2,max(maxx2(idx))+[.8 .8],'Color','k')
%         line(repmat((idx-2)*.22+2,2,1),max(maxx2(idx))+[.5 .5;.8 .8],'Color','k')
%         text((mean(idx)-2)*.22+2,max(maxx2(idx))+1.1,'*','Color','k','HorizontalAlignment','center')
%     end
% end
% 
% if p(2)<.05
%     maxx3 = max([maxx; maxx2],[],2);
%       line([.75  1.75; 1.25 2.25],max(maxx3)+[5 5;5 5],'Color','k');
%       
%     line([1 1],max(maxx3)+[5 7],'Color','k');
%     line([2 2],max(maxx3)+[5 7],'Color','k');
%     line([1 2],max(maxx3)+[7 7],'Color','k');
%         text(1.5,max(maxx3)+8,'*','Color','k','HorizontalAlignment','center')
% 
% end
%     
%     tmp = get(mm,'Position');
%     set(mm,...
%         'DefaultAxesFontSize',5,...
%         'DefaultLineLineWidth',2,...
%         'PaperUnits','points',...
%         'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
%         'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
%     saveas(mm,'FeatureResults/figs/Accuracy2WayANOVA.pdf','pdf');
%     
%     
%     
%     
%     
%     
% %Plot bitrates
% [p,atab,stats,terms] = anovan(data2wayB,{factor12wayB factor22wayB},'varnames',...
%     {'Type','Group'},'model','interaction')%,'display','off')
% [c1,m,h,nms] = multcompare(stats,'dimension',1,'display','off');
% [c2,m,h,nms] = multcompare(stats,'dimension',2,'display','off');
% 
%    
% mm = figure('Position',[100 100 300 400],'Name','P300 Bitrate across groups and tasks');
% barweb_AG([bit_mean_pat; bit_mean_con],[bit_err_pat; bit_err_con],.9,{'Patients','Controls'},[],[],[],[.3 .3 .3;.6 .6 .6;.9 .9 .9]);    hold on;
% plot([0 3],max([bit_err_pat+bit_mean_pat bit_err_con+bit_mean_con]),'--k'); 
% ylim([0 1.4*max([bit_err_pat+bit_mean_pat bit_err_con+bit_mean_con])]); xlim([.5 2.5]);
% legend('FS','RS','CS','Orientation','Horizontal')
% ylabel('Bitrate (bits/min)');
% 
% 
% if p(1)<.05 %significant effect of P300 type
%     maxx = bit_mean_pat+bit_err_pat;
%     for tt=find(sign(c1(:,3).*c1(:,5))==1)'
%         idx = c1(tt,1:2)
%         line((idx-2)*.22+1,max(maxx(idx))+[.8 .8],'Color','k')
%         line(repmat((idx-2)*.22+1,2,1),max(maxx(idx))+[.5 .5;.8 .8],'Color','k') 
%         text((mean(idx)-2)*.22+1,max(maxx(idx))+1.1,'*','Color','k','HorizontalAlignment','center')
%     end
%     maxx2 = bit_mean_con+bit_err_con;
%     for tt=find(sign(c1(:,3).*c1(:,5))==1)'
%         idx = c1(tt,1:2)
%         line((idx-2)*.22+2,max(maxx2(idx))+[.8 .8],'Color','k')
%         line(repmat((idx-2)*.22+2,2,1),max(maxx2(idx))+[.5 .5;.8 .8],'Color','k') 
%         text((mean(idx)-2)*.22+2,max(maxx2(idx))+1.1,'*','Color','k','HorizontalAlignment','center')
%     end
% end
% 
% if p(2)<.05
%     maxx3 = max([maxx; maxx2],[],2)
%       line([.75  1.75; 1.25 2.25],max(maxx3)+[2 2;2 2],'Color','k');
%       
%     line([1 1],max(maxx3)+[2 2.3],'Color','k');
%     line([2 2],max(maxx3)+[2 2.3],'Color','k');
%     line([1 2],max(maxx3)+[2.3 2.3],'Color','k');
%         text(1.5,max(maxx3)+2.6,'*','Color','k','HorizontalAlignment','center')
% 
% end
%     
%     tmp = get(mm,'Position');
%     set(mm,...
%         'DefaultAxesFontSize',5,...
%         'DefaultLineLineWidth',2,...
%         'PaperUnits','points',...
%         'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
%         'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
%     saveas(mm,'FeatureResults/figs/Bitrate2WayANOVA.pdf','pdf');
%     
%     

%%
save('OnlineResults.mat','P300_r','nms');




















