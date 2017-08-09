%Plots the average traces for each individual under different
%conditions.  Also gives the amplitude bar and the number of trials
%that went into each average.

%Variables input should be grouped into pairs of two, (P300, non-P300)
%(Left/Right) -- this way the counts can be displayed together

function [] = ALS_plottrace(pltname,vars,counts,plot_cts,ppl,pplrng,colors,...
    plot_sep,sepfactor,chans,channames,drange,sep_pad,dostat,comp)

text_pad = range(drange)/50;

tt = figure('Name',[pltname ' in ' channames{:}],'Position',[100 50 1200 700]);
kk=1;
for ii = 1:length(pplrng)
    for jj = 1:length(chans)
        if strcmp(ppl(ii),'P38') && sum(chans(jj)==[13 17 18 19])==1
            %Dont plot if P38 and one of these channels (they werent
            %placed in the correct location)
        else
            pp(jj)=subplot(1, length(chans), jj);
            for vv = 1:length(vars)
                plot(drange,squeeze(vars{vv}(pplrng(ii),chans(jj),:))-sepfactor(jj)*kk,'Color',colors{vv}); hold on;
            end
            set(gca,'Yticklabel',[]);
            if jj == 1
                text(drange(1)-.05,vars{vv}(pplrng(ii),chans(jj),1)-sepfactor(jj)*kk,...
                ppl{pplrng(ii)},'HorizontalAlignment','right');
            end
            if plot_cts == 1
                if jj == 1 && length(vars)>1
                    for vv = 1:2
                        text(drange(end)+text_pad,-sepfactor(jj)*kk-sepfactor(jj)*(vv-1.5)/2,...
                            ['n=' num2str(counts{vv}(pplrng(ii)))],...
                            'Color',colors{vv},'FontSize',8);
                    end
                end
                if jj == 2 && length(vars)>3
                    for vv = 1:2
                        text(drange(end)+text_pad,-sepfactor(jj)*kk-sepfactor(jj)*(vv-1.5)/2,...
                            ['n=' num2str(counts{vv+2}(pplrng(ii)))],...
                            'Color',colors{vv+2},'FontSize',8);
                    end
                end
                if jj == 3 && length(vars)>5
                    for vv = 1:2
                        text(drange(end)+text_pad,-sepfactor(jj)*kk-sepfactor(jj)*(vv-1.5)/2,...
                            ['n=' num2str(counts{vv+4}(pplrng(ii)))],...
                            'Color',colors{vv+4},'FontSize',8);
                    end
                end
                if jj == 4 && length(vars)>7
                    for vv = 1:2
                        text(drange(end)+text_pad,-sepfactor(jj)*kk-sepfactor(jj)*(vv-1.5)/2,...
                            ['n=' num2str(counts{vv+6}(pplrng(ii)))],...
                            'Color',colors{vv+6},'FontSize',8);
                    end
                end
            end
        end
    end
    kk=kk+1;
end
kk = sep_pad;
%Then plot grand averages
for jj = 1:length(chans)
    pp(jj)=subplot(1, length(chans), jj);
    for vv = 1:length(vars)
        if sum(chans(jj)==[13 17 18 19])==1
            %Average all participants except P38
            plot(drange,squeeze(nanmean(vars{vv}(~strcmp(ppl,'P38'),chans(jj),:)))-sepfactor(jj)*kk,'Color',colors{vv}); hold on;
        else
            plot(drange,squeeze(nanmean(vars{vv}(:,chans(jj),:)))-sepfactor(jj)*kk,'Color',colors{vv}); hold on;
        end
    end
    title(channames{jj});
    if jj==1
        text(drange(1)-.05,nanmean(vars{vv}(:,chans(jj),1))-sepfactor(jj)*kk,'Avg','HorizontalAlignment','right');
    end
    
    if plot_sep==1
    line([drange(end)+text_pad drange(end)+text_pad],[-sepfactor(jj)*(kk-.5) -sepfactor(jj)*(kk+.5)],...
        'Color','k','LineWidth',2,'clipping','off')
    text(drange(end)+2*text_pad, -sepfactor(jj)*(kk), sprintf(num2str(sepfactor(jj))), 'clipping','off');
    text(drange(end)+2*text_pad, -sepfactor(jj)*(kk)-sepfactor(jj)/2, '\muV', 'clipping','off');
    end
    ylim([-sepfactor(jj)*(kk+.5)-sepfactor(jj)/2 sepfactor(jj)/2]);
end


if dostat == 1
    for jj = 1:length(chans)
        pp(jj) = subplot(1,length(chans),jj);
        if sum(chans(jj)==[13 17 18 19])==1
            aa = squeeze(vars{comp(1)}(~strcmp(ppl,'P38'),chans(jj),:));
            aaW = squeeze(vars{comp(2)}(~strcmp(ppl,'P38'),chans(jj),:));
        else
            aa = squeeze(vars{comp(1)}(:,chans(jj),:));
            aaW = squeeze(vars{comp(2)}(:,chans(jj),:));
        end
        
        fs = 256;
        LPv = 20; %Filter at 20 hz, then downsample
        [bb2 aa2] = butter(5,LPv/(fs/2),'low');
        ad = filtfilt(bb2,aa2,aa');
        ad2 = downsample(ad,fix(fs/(LPv)))';
        awd = filtfilt(bb2,aa2,aaW');
        awd2 = downsample(awd,fix(fs/(LPv)))';
        ddrange = downsample(drange,fix(fs/LPv));
        clear at ar
        
        for iii = 1:size(ad2,2)
            [~,at(iii)] = ttest(ad2(:,iii),awd2(:,iii));
            ar(iii) = ranksum(ad2(:,iii),awd2(:,iii));
        end
        %Using the lillietest indicates that we can assume normality.
        w_sig=find(at<.05/size(ad2,2));
        
        for iii = 1:length(w_sig)
            line([ddrange(w_sig(iii))-mean(diff(ddrange))/2 ddrange(w_sig(iii))+mean(diff(ddrange))/2],...
                [-sepfactor(jj)*(kk+.5) -sepfactor(jj)*(kk+.5)],'Color','k','LineWidth',2);
        end
    end
end


linkaxes(pp);  xlim([drange(1) drange(end)]);
subplot(1,length(chans),2); xlabel('Time (seconds after stimulus)');


tmp = get(tt,'Position');
set(tt,...
    'DefaultAxesFontSize',5,...
    'DefaultLineLineWidth',2,...
    'PaperUnits','points',...
    'PaperSize',[tmp(3) tmp(4)],... was [11.5 8]
    'PaperPosition',[0 0 tmp(3) tmp(4)]); %was [0 0 11.5 8]
saveas(tt,['figures/' pltname '.pdf'],'pdf');

end