function [Data] = HeartbeatFilter(Data,State,Var,DoHB)
%Typically the heartbeat shows up in PO7, PO8, and Oz, but it can
%corrupt all channels.  Find the heartbeat in these three channels,
%by filtering out low frequency data, performing peak detection, creating a
%matched filter, refining the heartbeat selections, and subtracting the
%heartbeat template.

if DoHB == 1
    [bb,aa] = butter(4,[4 60]/(Var.fs/2));    
    HBdata = filtfilt(bb,aa,Data')';
    temp = mean(HBdata(6:8,:),1);
    
    
    %Expect roughly one beat per second
    expectedbeats = length(temp)/(Var.fs);
    i = 8;
    maxt = [];
    while size(maxt,1)<expectedbeats
        if i == 3
            disp('no heartbeat detected')
            break
        end
        [maxt, mint] = peakdet(temp,i*std(temp));
        i=i-1;
    end
    
    
    
    if ~isempty(maxt)
        maxt(maxt(:,1)<Var.fs|maxt(:,1)>size(Data,2)-Var.fs,:) = [];
        bpm = Var.fs*60./diff(maxt(:,1));
        %remove heartbeats greater than 200 bpm
        bpmi = bpm>200 | circshift(bpm>200,1);
        keepbeep = true(length(bpmi),1);
        counter = []; start = 0;
        for jj = 1:length(bpmi)
            if bpmi(jj) == 0
                if start == 1
                [xi,yi] = sort(counter(:,2));
                keepbeep(counter(yi(1:end-1),1)) = 0;
                counter = []; start = 0;
                end
            else
                counter = [counter; [jj maxt(jj,2)]];
                start = 1;
            end
        end
        if ~isempty(counter)
               [xi,yi] = sort(counter(:,2));
                keepbeep(counter(yi(1:end-1),1)) = 0;
                counter = []; start = 0;
        end
        maxt = maxt(find(keepbeep),:)
        bpm2 = Var.fs*60./diff(maxt(:,1));
        statschan = [mean(bpm2) std(bpm2)];
        maxval = maxt(:,2);
        maxloc = maxt(:,1);
        if Var.Figures_On==1
            figure
            yy(1) = subplot(3,4,1:3); plot(temp/max(temp),'Color',[.6 .6 .6]); hold on;
            scatter(maxloc, maxval/max(temp),100,'.k');
            subplot(3,4,4); hist(bpm2,200);
        end
    else
        statschan = [0 100];
    end
    
    if statschan(2)<40
        disp(['Preliminary heartbeat at ' num2str(statschan(1)) ' bpm, creating matched filter'])
        keepbeep2 = bpm2>statschan(2)-20 & bpm2<statschan(1)+20;
        tempind = repmat(maxloc(keepbeep2,1),1,Var.fs)+repmat(-(Var.fs/2-1):Var.fs/2,sum(keepbeep2),1);
        template = nanmean(temp(tempind));
        
        %         temp2 = conv(temp,template,'same');
        temp3 = conv(temp,fliplr(template),'same');
        
        %         figure
        %         xx(1) = subplot(311); plot(temp);
        %         xx(2) = subplot(312); plot(temp2);
        %         xx(3) = subplot(313); plot(temp3);
        %         linkaxes(xx,'x');
        
        i = 6;
        maxt3 = [];
        while size(maxt3,1)<expectedbeats
            [maxt3, mint3] = peakdet(temp3,i*std(temp3));
            i=i-1;
            if i == 0
                
                break
            end
        end
        if ~isempty(maxt3)
            maxt3(maxt3(:,1)<Var.fs|maxt3(:,1)>size(Data,2)-Var.fs,:) = [];
            bpm3 = Var.fs*60./diff(maxt3(:,1));
            %remove heartbeats greater than 200 bpm
            bpm3i = bpm3>200 | circshift(bpm3>200,1)
            keepbeep3 = true(length(bpm3i),1);
            counter = []; start = 0;
            for jj = 1:length(bpm3i)
                if bpm3i(jj) == 0
                    if start == 1
                        [xi,yi] = sort(counter(:,2));
                        keepbeep3(counter(yi(1:end-1),1)) = 0;
                        counter = []; start = 0;
                    end
                else
                    counter = [counter; [jj maxt3(jj,2)]];
                    start = 1;
                end
            end
            if ~isempty(counter)
                [xi,yi] = sort(counter(:,2));
                keepbeep3(counter(yi(1:end-1),1)) = 0;
                counter = []; start = 0;
            end
            maxt3 = maxt3(find(keepbeep3),:)
            bpm4 = Var.fs*60./diff(maxt3(:,1));
            statschan3 = [mean(bpm4) std(bpm4)];
            maxval3 = maxt3(:,2);
            maxloc3 = maxt3(:,1);
            if Var.Figures_On==1
                yy(2) = subplot(3,4,5:7); hold on; plot(temp3/max(temp3),'Color',[1 .5 .5]);
                scatter(maxloc3, temp3(maxloc3)/max(temp3),100,'.r');
                subplot(3,4,8); plot(template,'k');
            end
            
            Data2 = Data;
            tempind = repmat(maxt3(:,1),1,Var.fs)+repmat(ceil(-Var.fs/2+1:Var.fs/2),size(maxt3,1),1);
            wnd = window(@hanning,fs);
            for i = 1:8
                tempd = reshape(HBdata(i,tempind'),fs,size(tempind,1));
                %remove artifactual epochs
                bade = var(tempd)>5*abs(mean(var(tempd)));
                tempd(:,bade) = [];
                signature(:,i) = wnd.*mean(tempd,2);
                
                for j = 1:size(maxt3,1)
                    Data2(i,maxt3(j,1)+ceil(-Var.fs/2+1:Var.fs/2)) = Data2(i,maxt3(j,1)+ceil(-Var.fs/2+1:Var.fs/2))-signature(:,i)';
                end
            end
            
            if Figures_On==1
                subplot(3,4,12); plot(signature);
                yy(3) = subplot(3,4,9:11);  plot(mean(Data(6:8,:),1),'k'); hold on; 
                plot(mean(Data2(6:8,:),1),'Color',[.5 .5 .5]);
                plot(mean(Data(1,:),1)-50,'b'); hold on; plot(mean(Data2(1,:),1)-50,'Color',[.5 .5 1]);
                linkaxes(yy,'x');
            end
            
        else
            disp('no heartbeat filter created');
        end
        
        
        
        %         %I dont think ICA is ideal for this situation.  First, it takes
        %         %long.  Second there are only eight channels.
        %         [weights, sphere] = runica(Data);
        %         temp = weights*Data;
        %         figure
        %         plot(temp');
        %
        %         for i = 1:8
        %             [maxt, mint] = peakdet(temp(i,:),5*std(temp(i,:)));
        %             bpm = fs*60./diff(maxt(:,1));
        %             %remove heartbeats greater than 100 bpm
        %             keepbeep = bpm<=100;
        %             bpm(~keepbeep,:) = [];
        %             scatter(maxt(keepbeep,1),maxt(keepbeep,2));
        %             subplot(1,4,4); hist(bpm,200);
        %             statsica(i,:) = [mean(bpm) std(bpm)];
        %         end
        %         [~, candidate] = min(statsica(:,2));
        %         if abs((statsica(candidate,1)-statschan(1))/statschan(1)*100)<5 %bpm
        %             %from the candidate ica channel is less than 5% off the bpm in the
        %             %PO/O channels.  Remove this signal from data.
        %             invweights = inv(weights);
        %             invweights(:,candidate) = 0;
        %             temp2 = invweights*temp;
        
        
        if Var.Figures_On == 1
            figure
            xx(1) = subplot(211); plot(Data'); hold on; plot(State.StimulusCodeUNArt,'k','LineWidth',1.5);
            xx(2) = subplot(212); plot(Data2'); hold on; plot(State.StimulusCode,'k','LineWidth',1.5);
            linkaxes(xx);
        end
        Data = Data2;
    else
        disp('Heartbeat signal too noisy... skipping filter');
    end
else
    disp('Heartbeat filter disabled')
end


clear HBdata
