function gazetracking(ppl)
close all

REF = 0;
fltext = '';
drange = round(.45*250):round(1.2*250);
pplrng_t = 1:length(ppl);

aa = 1;
cfm = [];



gf_acc_x_tot = cell(length(ppl),1); %32 trials over two sessions
gf_acc_y_tot = cell(length(ppl),1);
gf_var_x_tot = cell(length(ppl),1);
gf_var_y_tot = cell(length(ppl),1);
gf_inv_tot = cell(length(ppl),1);
sessctr = cell(length(ppl),1);
pplctr = cell(length(ppl),1);

curppl = [];
for ii = pplrng_t
    
    try
        load(['results\'...
            ppl{ii} fltext '_P300feats_Ref' num2str(REF)]);
        
        xx = figure('Name',ppl{ii},...
            'Position',[100 100 1200 500],'visible','on');
        nogo = 0;
    catch
        disp(['No data for ' ppl{ii}]);
        gf_acc_x_tot{ii} = NaN(1,32);
        gf_acc_y_tot{ii} = NaN(1,32);
        gf_var_x_tot{ii} = NaN(1,32);
        gf_var_y_tot{ii} = NaN(1,32);
        gf_inv_tot{ii} = NaN(1,32);
        nogo = 1;
    end
    
    if ii == 5
       hi = 1 
    end
    
    %Convert all values from pixels to cm.  Screen used
    %was 1920x1080.  Dimensions of monitor are 47.8cm by 26.9cm.
    xpixconv = 1920/47.7;
    ypixconv = 1080/26.9;
    
    
    
    if ~nogo
        %             gf_acc_x_temp = NaN(1,32); gf_acc_y_temp = NaN(1,32);
        %             gf_var_x_temp = NaN(1,32); gf_var_y_temp = NaN(1,32);
        %             gf_inv_temp = NaN(1,32);
        for sess = 1:size(GAZE_FULL,1)
            iCol = mod(sess-1,4)+1
            iRow = ceil(sess/4)
            SUPERPLOT(2,4,iRow, iCol); axis off;
            %                 rectangle('Position',[1600 -1900 800 800]./...
            %                     [xpixconv ypixconv xpixconv ypixconv],'EdgeColor',[.8 .8 .8])
            %                 line([1600 2400]./(xpixconv*[1 1]),[-1220 -1220]./(ypixconv*[1 1]),'Color',[.8 .8 .8])
            %
            
         
            
            gf = GAZE_FULL(sess,:)';
            gf = gf(cellfun(@(x) ~isempty(x),gf));
            if ~isempty(gf)
                gf_targs = []; gf_acc = []; gf_var = []; gf_inv = [];
                for rnn = 1:length(gf)
                    temp_gf = gf{rnn};
                    if isempty(temp_gf)
                        gf_targs = [gf_targs; NaN(4,4)];
                        gf_acc = [gf_acc; NaN(4,2)];
                        gf_var = [gf_var; NaN(4,2)];
                        gf_inv = [gf_inv; 100*ones(4,1)];
                    else
                        gf_targs = [gf_targs; gf{rnn}.targetbox./...
                            repmat([xpixconv ypixconv xpixconv ypixconv],...
                            size(gf{rnn}.targetbox,1),1)];
                        gf_acc = [gf_acc; gf{rnn}.acc./...
                            repmat([xpixconv ypixconv],...
                            size(gf{rnn}.targetbox,1),1)];
                        gf_var = [gf_var; gf{rnn}.var./...
                            (repmat([xpixconv ypixconv],...
                            size(gf{rnn}.targetbox,1),1).^2)];
                        gf_inv = [gf_inv; gf{rnn}.invalid];
                    end
                end
                
                gf_acc_x_tot{ii,sess} = gf_acc(:,1)';
                gf_acc_y_tot{ii,sess} = gf_acc(:,2)';
                gf_var_x_tot{ii,sess} =  gf_var(:,1)';
                gf_var_y_tot{ii,sess} =  gf_var(:,2)';
                sessctr{ii,sess} = repmat(sess,1,size(gf_acc,1));
                pplctr{ii,sess} = repmat(ii,1,size(gf_acc,1));
                tmp = gf_acc(:,1);
                tmp(~isnan(gf_acc(:,1)))=1;
                gf_inv_tot{ii,sess} = tmp'.*gf_inv';
                
                
                
                temp = sqrt(gf_acc_x_tot{ii,sess}.^2+gf_acc_y_tot{ii,sess}.^2)
                AvgEye.acc{ii,sess} = nanmean(temp);
                temp = sqrt(gf_var_x_tot{ii,sess}.^2+gf_var_y_tot{ii,sess}.^2);
                AvgEye.pre{ii,sess} = nanmean(temp);
                AvgEye.inv{ii,sess} = nanmean(gf_inv_tot{ii,sess});
                AvgEye.num{ii,sess} = sum(~isnan(gf_inv_tot{ii,sess}),2);
                
                
                
                
                
                [A, ~, IA] = unique(gf_targs,'rows');
                for ia = 1:max(IA)
                    numa(ia) = sum(IA==ia);
                end
                
                
                hold on;
                
                ctr = [(gf_targs(:,1)+gf_targs(:,3))/2,...
                    (gf_targs(:,2)+gf_targs(:,4))/2];
                
                
                [um,ix,jx] = unique(gf_targs,'rows');
                dupeinds = setdiff([1:size(gf_targs,1)],ix);
                for gt = 1:size(gf_targs,1)
                    pts = [gf_targs(gt,1),gf_targs(gt,1),gf_targs(gt,3),gf_targs(gt,3);...
                        gf_targs(gt,2),gf_targs(gt,4),gf_targs(gt,4),gf_targs(gt,2)];
                    %Made the fill only half solid if there was 100% invalid.
                    %This way, you can see if there were multiple trials which
                    %were invalid.
                    patch(pts(1,:),pts(2,:),'r','FaceAlpha',gf_inv(gt)/200);
                    
                    
                    ellipse(sqrt(gf_var(gt,1)),sqrt(gf_var(gt,2)),[],ctr(gt,1)+gf_acc(gt,1),...
                        ctr(gt,2)+gf_acc(gt,2),'k');
                    if ismember(gt,dupeinds)
                        text(pts(1,2)+.1,pts(2,2)+.3,num2str(gt))
                    else
                        text(pts(1,1)+.1,pts(2,1)-.3,num2str(gt))
                    end
                end
                quiver(ctr(:,1),ctr(:,2),gf_acc(:,1),gf_acc(:,2),0,'Color','k','LineWidth',1)
                axis equal;
                if mean(ctr(:,1))<75
                    ylim([-1940 -1320]/ypixconv);
                    xlim([2390 3210]/xpixconv);
                    rectangle('Position',[2400 -1925.2 800 595.2]./...
                        [xpixconv ypixconv xpixconv ypixconv],'EdgeColor',[.8 .8 .8])
                    line([2400 3200]./(xpixconv*[1 1]),[-1419 -1419]./(ypixconv*[1 1]),'Color',[.8 .8 .8])
                    text(70,-34,['Session' num2str(sess)],'FontSize',18,'HorizontalAlignment','center');
                    
                else
                    ylim([-1940 -1320]/ypixconv);
                    xlim([2890 3710]/xpixconv);
                    rectangle('Position',[2900 -1925.2 800 595.2]./...
                        [xpixconv ypixconv xpixconv ypixconv],'EdgeColor',[.8 .8 .8])
                    line([2900 3700]./(xpixconv*[1 1]),[-1419 -1419]./(ypixconv*[1 1]),'Color',[.8 .8 .8])
                    text(80,-34,['Session' num2str(sess)],'FontSize',18,'HorizontalAlignment','center');
                   
                end
                
                
            else
                SUPERPLOT(2,4,iRow, iCol); axis off;
                text(.2,.8,'No Eye Gaze Data');
                AvgEye.acc{ii,sess} = NaN;
                AvgEye.pre{ii,sess} = NaN;
                AvgEye.inv{ii,sess} = NaN;
                AvgEye.num{ii,sess} = NaN;
            end
        end
        
        
        tmp = get(xx,'Position');
        set(xx,...
            'DefaultAxesFontSize',5,...
            'DefaultLineLineWidth',2,...
            'PaperUnits','points',...
            'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
            'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
        saveas(xx,['figures/' ppl{ii} '_Gaze_new.pdf'],'pdf');
        close(xx);
        
        curppl = [curppl ppl(ii)];
    end
    [~,AvgEye.minsessacc(ii)] = min(cell2mat(AvgEye.acc(ii,:)),[],2);
    %For accuracy, precision, and invadility, only calculate these states from
    %the session with the lowest average error.
    
    AvgEye.minacc(ii) = AvgEye.acc{ii,AvgEye.minsessacc(ii)};
    AvgEye.minpre(ii) = AvgEye.pre{ii,AvgEye.minsessacc(ii)};
    AvgEye.mininv(ii) = AvgEye.inv{ii,AvgEye.minsessacc(ii)};
    AvgEye.minnum(ii) = AvgEye.num{ii,AvgEye.minsessacc(ii)};
end







all_gd = {gf_acc_x_tot',gf_acc_y_tot',...
    gf_var_x_tot', gf_var_y_tot',...
    gf_inv_tot', sessctr', pplctr'};
all_gd = cellfun(@(x) x(:),all_gd,'UniformOutput',false);

clear dtt
for j = 1:length(all_gd)
    tmpgd = all_gd{j};
    maxtrial = max(cellfun(@length,tmpgd));
    dtt{j} = NaN(length(tmpgd),maxtrial);
    for i= 1:size(dtt{j},1)
        dtt{j}(i,1:length(tmpgd{i})) = tmpgd{i};
    end
end

% 
% all_gd = cellfun(@(x) x(:),all_gd,'UniformOutput',false)
% all_gd_buff = cellfun(@(x) [x x(:,end)],all_gd,'UniformOutput',false);
% all_gd_buff = cellfun(@(x) [x; x(end,:)],all_gd_buff,'UniformOutput',false);
all_gd_label = {'Horizontal Gaze Error (cm)',...
    'Vertical Gaze Error (cm)',...
    'Horizontal Gaze Variance (cm^2)',...
    'Vertical Gaze Variance (cm^2)'};
%Limit the range to 80 pixels, this is the height and width of the
%targetbox
all_gd_rng = {[-5 5],[-5 5],[0 5],[0 5]};
mm = figure('Name','GazeStats','Position',[100 100 1200 500]);
for jj = 1:4
    
    subplot(1,4,jj);
    h = pcolor(dtt{jj}); shading flat; cb = colorbar; caxis(all_gd_rng{jj})
    ylabel(cb,all_gd_label{jj});
    axis off;
    for i =1:size(gf_acc_x_tot,1)
        line([-3 -3],(i-1)*10+[1 10]+.5,'Color','k','Clipping','off',...
            'LineWidth',2)
        text(-4,(i-1)*10+5.5,ppl(i),'HorizontalAlignment','right');
        text(zeros(2,1),(i-1)*10+[1 10]+.5,{'1','10'},...
            'HorizontalAlignment','right')
    end
end


tmp = get(mm,'Position');
set(mm,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(mm,['figures/GazeStats.pdf'],'pdf');
close(mm);




for ii = pplrng_t
    mm=figure; hold on;
    cmp = parula(10);
    for jj = 1:length(gf_var_x_tot) %sessions
        ctr = [(gf_targs(1,1)+gf_targs(1,3))/2,...
            (gf_targs(1,2)+gf_targs(1,4))/2];
        
        pts = [gf_targs(1,1),gf_targs(1,1),gf_targs(1,3),gf_targs(1,3);...
            gf_targs(1,2),gf_targs(1,4),gf_targs(1,4),gf_targs(1,2)];
        %Made the fill only half solid if there was 100% invalid.
        %This way, you can see if there were multiple trials which
        %were invalid.
        patch(pts(1,:),pts(2,:),'r','FaceAlpha',gf_inv(gt)/200,'LineWidth',2);
        
        
        h = ellipse(sqrt(nanmean(gf_var_x_tot{ii,jj})),sqrt(nanmean(gf_var_y_tot{ii,jj})),...
            [],ctr(1)+nanmean(gf_acc_x_tot{ii,jj}),...
            ctr(2)+nanmean(gf_acc_y_tot{ii,jj}),cmp(jj,:));
        set(h,'LineWidth',2);
        
        text(ctr(1)+nanmean(gf_acc_x_tot{ii,jj})-.8*sqrt(nanmean(gf_var_x_tot{ii,jj})),...
            ctr(2)+nanmean(gf_acc_y_tot{ii,jj})+.8*sqrt(nanmean(gf_var_y_tot{ii,jj})),num2str(jj),...
            'HorizontalAlignment','center','Color',cmp(jj,:),'Fontsize',18)
        
        
        
        quiver(ctr(1),ctr(2),nanmean(gf_acc_x_tot{ii,jj}),...
            nanmean(gf_acc_y_tot{ii,jj}),0,'Color',cmp(jj,:),'LineWidth',2)
        
    end
    axis equal;
    xlabel('X position (cm)');
    ylabel('Y position (cm)');
    set(gca,'FontSize',18);
    tmp = get(mm,'Position');
    set(mm,...
        'DefaultAxesFontSize',5,...
        'DefaultLineLineWidth',2,...
        'PaperUnits','points',...
        'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
        'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
    saveas(mm,['figures/EyeAccBySession_' num2str(ii) '.pdf'],'pdf');
    close(mm);
end





for i = 1:10
    gf_acc_x_GA{i} = abs(cat(2,gf_acc_x_tot{[1 3 5],i}));
    gf_acc_y_GA{i} = abs(cat(2,gf_acc_y_tot{[1 3 5],i}));
    gf_var_x_GA{i} = abs(cat(2,gf_var_x_tot{[1 3 5],i}));
    gf_var_y_GA{i} = abs(cat(2,gf_var_y_tot{[1 3 5],i}));
end


errmag = sqrt(cellfun(@nanmean,gf_acc_x_GA).^2+cellfun(@nanmean,gf_acc_y_GA).^2);
hh = figure; hold on; k =1; axis equal;
cmp = parula(10);
for i = 0:324/9:324
    [x,y] = pol2cart(pi*i/180,errmag(k))
    quiver(0,0,x,y,'Color',cmp(k,:),'LineWidth',2)
    h = ellipse(nanmean(sqrt((gf_var_x_GA{k}))),nanmean(sqrt(gf_var_y_GA{k})),...
        [],x,y,cmp(jj,:));
    set(h,'LineWidth',2,'Color',cmp(k,:));
    k=k+1;
end
tmp = get(hh,'Position');
set(hh,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(hh,['figures/EyeAccBySession_' num2str(ii) '.pdf'],'pdf');
close(hh);
%
%
%
% pFullIndv = [pFullIndv {{AvgEye.acc{1}(prng_S2) AvgEye.pre{1}(prng_S2) AvgEye.inv{1}(prng_S2) ...
%     AvgEye.acc{2}(prng_S2) AvgEye.pre{2}(prng_S2) AvgEye.inv{2}(prng_S2) ...
%     AvgEye.acc{3}(prng_S2) AvgEye.pre{3}(prng_S2) AvgEye.inv{3}(prng_S2)}}];
% pFullIndv_L = [pFullIndv_L {'Eye Tracking'}];
% pFullIndv_Ls = [pFullIndv_Ls {{'Gaze_{FSerr}','Gaze_{FSvar}','Gaze_{FSinv}',...
%     'Gaze_{RSerr}','Gaze_{RSvar}','Gaze_{RSinv}',...
%     'Gaze_{CSerr}','Gaze_{CSvar}','Gaze_{CSinv}'}}];
%
%
% cFullIndv = [cFullIndv {{AvgEye.acc{1}(crng_S2) AvgEye.pre{1}(crng_S2) AvgEye.inv{1}(crng_S2) ...
%     AvgEye.acc{2}(crng_S2) AvgEye.pre{2}(crng_S2) AvgEye.inv{2}(crng_S2) ...
%     AvgEye.acc{3}(crng_S2) AvgEye.pre{3}(crng_S2) AvgEye.inv{3}(crng_S2)}}];
% cFullIndv_L = [cFullIndv_L {'Eye Tracking'}];
% cFullIndv_Ls = [cFullIndv_Ls {{'Gaze_{FSerr}','Gaze_{FSvar}','Gaze_{FSinv}',...
%     'Gaze_{RSerr}','Gaze_{RSvar}','Gaze_{RSinv}',...
%     'Gaze_{CSerr}','Gaze_{CSvar}','Gaze_{CSinv}'}}];
%
%
%
% %No differences for any of the tasks between patients and controls in gaze
% %accuracy, precision, or validity of data.
% EyeAccDiff = cellfun(@(x) ranksum(x(prng_S2),x(crng_S2)),AvgEye.acc);
% EyePreDiff = cellfun(@(x) ranksum(x(prng_S2),x(crng_S2)),AvgEye.pre);
% EyeInvDiff = cellfun(@(x) ranksum(x(prng_S2),x(crng_S2)),AvgEye.inv);
%
% %C02 and P66 are consistant outliers.  Remove and redo testing
% olrp = find(strcmp(ppl_S2,'P66'));
% olrc = find(strcmp(ppl_S2,'C02'));
% prng_S2o=prng_S2; prng_S2o(olrp)=0;
% crng_S2o=crng_S2; crng_S2o(olrc)=0;
% EyeAccDiff_withoutOutliers = cellfun(@(x) ranksum(x(prng_S2o),x(crng_S2o)),AvgEye.acc);
% EyePreDiff_withoutOutliers = cellfun(@(x) ranksum(x(prng_S2o),x(crng_S2o)),AvgEye.pre);
% EyeInvDiff_withoutOutliers = cellfun(@(x) ranksum(x(prng_S2o),x(crng_S2o)),AvgEye.inv);
%
%
% nRows = 3;
% nCols = 3;
% sMar = 0.2;
% tMar = 0.2;
% bMar = 0.2;
% cSpace = 0.02;
% vSpace = 0.05;
% height = (1-tMar-bMar-(nRows-1)*vSpace)/(nRows);
% width = (1-(nCols-1)*cSpace-2*sMar)/(nCols);
% %Show error and precision of patients and controls
% mm = figure('Name','Error and Precision of eye tracking in patients and controls','Position',[100 100 600 800]);
% kk=1;
% dta = {AvgEye.acc AvgEye.pre AvgEye.inv}; lbl1 = {'Error (cm)','Variance (cm^2)','Invalidity (% missing)'};
% sigg = [EyeAccDiff' EyePreDiff' EyeInvDiff'];
% for j = 1:3 %accuracy, precision, invalid
%     for i = 1:3 %Types of spellers
%
%         vB = (1-tMar)-(j)*(height+vSpace);
%         hL = sMar+(i-1)*(width+cSpace);
%         subplot('Position',[hL, vB, width, height]);
%         pt = dta{j}{i}(prng_S2);
%         ct = dta{j}{i}(crng_S2);
%
%         mtt = [nanmean(pt) nanmean(ct)];
%         ptt = [nanstd(pt)/sqrt(length(pt)) nanstd(ct)/sqrt(length(ct))];
%
%         if j==3 && i==2
%             boxplot([pt;ct],[ones(length(pt),1);2*ones(length(ct),1)],'colors',...
%                 [.3 .3 .3;.7 .7 .7],'jitter',.5,'labels',{'Patients','Controls'},'outliersize',4);
%         else
%             %       boxplot([pt;ct],[ones(length(pt),1);2*ones(length(ct),1)],'plotstyle','compact','colors',...
%             %           [.3 .3 .3;.7 .7 .7],'jitter',0,'labels',{'',''});
%             boxplot([pt;ct],[ones(length(pt),1);2*ones(length(ct),1)],'colors',...
%                 [.3 .3 .3;.7 .7 .7],'jitter',.5,'labels',{'',''},'outliersize',4);
%         end
%         %     barweb_AG(mtt,ptt,.9,[], [], [], [], summer);
%         %     xlim([.7 1.3])
%         %     if i==3 && j==3
%         %         legend('Patients','Controls');
%         %     end
%         %
%         %
%
%
%         %TO "REMOVE" two outliers (P66 and C02), set these limits.
%         if j==1
%             ylim([0 5]);
%         elseif j==2
%             ylim([-1 37]);
%         elseif j==3
%             ylim([-1 25]);
%         end
%         yext = range(get(gca,'ylim'));
%
%         if sigg(i,j)<.05
%             text(1.5,max([mtt+ptt])+.13*yext,'*','FontSize',18,'HorizontalAlignment','center');
%         end
%         if j==1
%             title(SpellerType{i})
%         end
%         if i==1
%             ylabel(lbl1{j},'FontWeight','bold');
%         else
%             set(gca,'Yticklabel',[]);
%         end
%         kk=kk+1;
%     end
% end
% tmp = get(mm,'Position');
% set(mm,...
%     'DefaultAxesFontSize',5,...
%     'DefaultLineLineWidth',2,...
%     'PaperUnits','points',...
%     'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
%     'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
% saveas(mm,['FeatureResults/figs/GroupTrackingDiff_withoutOutliers.pdf'],'pdf');
% close(mm);
%
%
%
% %Does the presence of glasses or ocular issues affect eye tracking
% EyeFac(:,1) = EyeHist.Gla(prng_S2R);
% EyeFac(:,2) = sum(EyeHist.Full(prng_S2R,1:3),2)>0;
% EyeFac_l = {'Glasses','Ocular Issues'};
% Occ_l = {'Tracking Error','Tracking Variance'};
% clear pv
% p_i = 1;
% mm = figure
% for ffc = 1:2 %Factors (glasses, eye issues)
%     for occ = 1:2 %outcome (error or variance of gaze)
%         subplot(2,2,p_i);
%         p_i = p_i+1;
%         for tpp= 1:3 %P300 type
%             if occ == 1
%                 temp = AvgEye.acc{tpp}(prng_S2);
%             elseif occ == 2
%                 temp = AvgEye.pre{tpp}(prng_S2);
%             end
%             tempm(ffc,occ,tpp,:) = [nanmean(temp(EyeFac(:,ffc)==0)) nanmean(temp(EyeFac(:,ffc)==1))];
%             temps(ffc,occ,tpp,:) = [nanstd(temp(EyeFac(:,ffc)==0))/sqrt(sum(EyeFac(:,ffc)==0))...
%                 nanstd(temp(EyeFac(:,ffc)==1))/sqrt(sum(EyeFac(:,ffc)==1))];
%             pv(ffc,occ,tpp) = ranksum(temp(EyeFac(:,ffc)==0),temp(EyeFac(:,ffc)==1))
%
%             barweb_AG(squeeze(tempm(ffc,occ,:,:))',squeeze(temps(ffc,occ,:,:))',...
%                 .9,{'Yes','No'},EyeFac_l{ffc},[],Occ_l{occ},[.3 .3 .3;.6 .6 .6; .9 .9 .9]);
%             if pv(ffc,occ,tpp)<.2
%                 hold on;
%                 disp('what')
%                 text(1.5,.2*tpp,'*','FontSize',18,'Color',tpp*[.3 .3 .3])
%             end
%         end
%
%     end
% end
% tmp = get(mm,'Position');
% set(mm,...
%     'DefaultAxesFontSize',5,...
%     'DefaultLineLineWidth',2,...
%     'PaperUnits','points',...
%     'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
%     'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
% saveas(mm,['FeatureResults/figs/EyeHistVStracking.pdf'],'pdf');
% close(mm);
%
%
%
% % prsn = 22; %P78
% % lmts = [-30 30; -30 30; -1 30; -1 30; -2 100];
% prsn = 12; % P56
% lmts = [-2 2; -2 2; -.01 .50; -.01 .50; -2 100];
% cmps = {[[repmat(linspace(0,1,50)',1,2) ones(50,1)];[ones(50,1) repmat(linspace(1,0,50)',1,2)]],...
%     [[repmat(linspace(0,1,50)',1,2) ones(50,1)];[ones(50,1) repmat(linspace(1,0,50)',1,2)]],...
%     [ones(100,1) repmat(linspace(1,0,100)',1,2)],...
%     [ones(100,1) repmat(linspace(1,0,100)',1,2)],...
%     [ones(100,1) repmat(linspace(1,0,100)',1,2)]};
% ylbl = {[{'X-error'} {'(cm)'}],[{'Y-error'} {'(cm)'}],...
%     [{'X-dev'} {'(cm^2)'}],[{'Y-dev'} {'(cm^2)'}],[{'Invalid'} {'Trials (%)'}]};
% mm=figure('Name',['Example data for ' ppl_S2{prsn}],'Position',[100 100 300 200]);
% for i = 1:5
%     if i == 5
%         tmpp = all_gd{1}{i}(prsn,1:16);
%         tmpp(isnan(tmpp))=100;
%         ax(i) = subplot(5,1,i); imagesc(tmpp,lmts(i,:));
%         xlabel('Trial Number','FontSize',10);
%     elseif i == 1
%         ax(i) = subplot(5,1,i); imagesc(all_gd{1}{i}(prsn,1:16),lmts(i,:));
%         text(1:16,.15*ones(16,1),{'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16'},...
%             'HorizontalAlignment','center','FontSize',8);
%     else
%         ax(i) = subplot(5,1,i); imagesc(all_gd{1}{i}(prsn,1:16),lmts(i,:));
%     end
%     colorbar; set(gca,'Xticklabel',[],'Yticklabel',[]);
%     text(.3,1.3,sprintf('%s\n',ylbl{i}{:}),'rotation',0','HorizontalAlignment','right','FontSize',10);
%     cmp = cmps{i};
%     colormap(ax(i),cmp)
%     cmp = colormap(ax(i));
%     cmp(1,:) = [.5 .5 .5];
%     colormap(ax(i),cmp);
% end
% tmp = get(mm,'Position');
% set(mm,...
%     'DefaultAxesFontSize',5,...
%     'DefaultLineLineWidth',2,...
%     'PaperUnits','points',...
%     'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
%     'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
% saveas(mm,['figures/ExampleGazeStats_' ppl_S2{prsn} '.pdf'],'pdf');
% close(mm);
