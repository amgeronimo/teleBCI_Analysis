function [Var] = ParseLogFile(Data,State,Var,InputText,TTSpell)




if Var.Speller == 1  %Use the myfile.txt associated with CB speller runs to
    %find stimulus codes and target associations
    
    for i = 1:length(InputText)
        RunEnd(i) = ~isempty(strfind(InputText{i},'******************'));
    end
    Runind = find(RunEnd==1);
    SCindF = []; SRindF = []; NEXTmindF = []; GoalLindF = [];
    Var.ChosenLetterF = []; Var.GoalLetterF = []; INITmindF = []; Var.GoalTrialF = 0;
    NEXTgoalletterF = [];
    for rr = 1:length(Runind);
        if rr == 1;
            RunRange = 1:Runind(1);
        else
            RunRange = Runind(rr-1):Runind(rr);
        end
        clear Begin SCind SRind INITmind NEXTmind ChosenLind GoalLind
        for ll = RunRange
            GoalLind(ll) = ~isempty(strfind(InputText{ll},'Goal Text,'));
            ChosenLind(ll) = ~isempty(strfind(InputText{ll},'Selected Text'));
            Begin(ll) = ~isempty(strfind(InputText{ll},'*****START'));
            SCind(ll) = ~isempty(strfind(InputText{ll},'SC'));
            SRind(ll) = ~isempty(strfind(InputText{ll},'SR'));
            INITmind(ll) = ~isempty(strfind(InputText{ll},'INIT_mSequence'));
            NEXTmind(ll) = ~isempty(strfind(InputText{ll},'NEXT_mSequence'));
        end
        %Ignore the last SC, SR, NEXTmind, and 2 GoalLetters
        GoalLind2 = find(GoalLind);
        ChosenLind2 = find(ChosenLind)+1; %this is the distance between the selected letter and the repeat of the next goal letter.
        trialkeep = ~ismember(GoalLind2,ChosenLind2);
        trialkeep(end)=0;
        GoalLindF = [GoalLindF (GoalLind2(trialkeep))];
        GoalLetter = InputText(GoalLind2(trialkeep));
        
        [GoalTrialbeg GoalTrialend] = cellfun(@(x) regexp(x,'\s\d*\s='),GoalLetter,'UniformOutput',false);
        GoalTrial = cellfun(@(x,y,z) str2num(x(y+1:z-2)),GoalLetter,GoalTrialbeg,GoalTrialend);
        stind = regexp(GoalLetter,'('); endind = regexp(GoalLetter,')');
        for xy = 1:length(stind)
            GoalLetter{xy} = GoalLetter{xy}(stind{xy}+1:endind{xy}-1);
        end
        
        
        GoalTrial = GoalTrial+max(Var.GoalTrialF);
        Var.GoalTrialF = [Var.GoalTrialF; GoalTrial];%This
        %vector is used to reconcile the number of trials and the number of letters spelled
        %             GoalLetter = TTSpell_cell{rr}(GoalTrial)';
        Var.GoalLetterF = [Var.GoalLetterF; GoalLetter];
        ChosenLetter = InputText(ChosenLind(1:end-1));
        ChosenLetter = cellfun(@(x) x(18:end),ChosenLetter,'UniformOutput',false);
        
        %             ChosenLetter = ChosenLetter(:,18);
        Var.ChosenLetterF = [Var.ChosenLetterF; ChosenLetter];
        Begind = find(Begin==1);
        SCind = find(SCind==1)';
        FSCi = find(SCind < Begind);
        SCindF = [SCindF; SCind(FSCi(end):end-1)];
        SRind = find(SRind==1)';
        FSRi = find(SRind < Begind);
        SRindF = [SRindF; SRind(FSRi(end):end-1)];
        INITmind = find(INITmind==1)';
        INITmindF = [INITmindF; INITmind(1:end-1)];
        NEXTmind = find(NEXTmind==1)';
        omitN = Var.NumSequences{rr}:Var.NumSequences{rr}:Var.NumTrials{rr}*Var.NumSequences{rr};
        NEXTmind(omitN) = [];
        NEXTmindF = [NEXTmindF; NEXTmind];
        tmpind = repmat(GoalTrial,1,Var.NumSequences{rr})';
        tmpind(omitN) = [];
        NEXTgoalletterF = [NEXTgoalletterF; tmpind(:)];
    end
    Var.GoalTrialF = Var.GoalTrialF(2:end);
    
    
    
    
    %         for j = 1:NumSequences-1
    %             tempNEXTm = InputText{NEXTmindF((NumSequences-1)*(i-1)+j)}(16:end);
    %             tempNEXTm = regexp(tempNEXTm,' ','split')'; tempNEXTm(cellfun(@isempty,tempNEXTm)) = [];
    %             NEXTm{j,i} = tempNEXTm;
    %         end
    %     end
    %     if isempty(NEXTm)
    %         mSeq = INITm;
    %     else
    %         mSeq = cat(1,INITm,NEXTm);
    %     end
    
    
    
    
    NEXTm = cell(1,length(Var.GoalLetterF));
    for i = 1:length(Var.GoalLetterF)
        tempSC = InputText{SCindF(i)}(4:end);
        tempSC = regexp(tempSC,' ','split')'; tempSC(cellfun(@isempty,tempSC)) = [];
        Var.SC{i} = tempSC;
        tempSR = InputText{SRindF(i)}(4:end);
        tempSR = regexp(tempSR,' ','split')'; tempSR(cellfun(@isempty,tempSR)) = [];
        Var.SR{i} = tempSR;
        tempINITm = InputText{INITmindF(i)}(16:end);
        tempINITm = regexp(tempINITm,' ','split')';
        tempINITm = cellfun(@(x) str2num(x),(tempINITm(1:end-1,:)));
        INITm{i} = tempINITm;
        
        %tempNEXTm = InputText{NEXTmindF((NumSequences{rr}-1)*(i-1)+j)}(16:end);
        tempNEXTm = InputText(NEXTmindF(NEXTgoalletterF==i));
        tempNEXTm = cellfun(@(x) x(16:end), tempNEXTm,'UniformOutput',false);
        tempNEXTm = regexp(tempNEXTm,'\s','split')';
        tempNEXTm = cat(1,tempNEXTm{:})';
        if ~isempty(tempNEXTm)
            tempNEXTm = cellfun(@(x) str2num(x),(tempNEXTm(1:end-1,:)));
        end
        NEXTm{i} = tempNEXTm;
        
        
        if isempty(NEXTm{i})
            mSeq{i} = INITm{i};
        else
            mSeq{i} = cat(2,INITm{i},NEXTm{i});
        end
    end
    %%%%%%%%%%%
    
elseif Var.Speller == 2
    
    for i = 1:sum(cell2mat(Var.NumTrials))
        Var.SC{i} = 1:Var.NC;
        Var.SR{i} = (Var.NC+1):2*Var.NC;
    end
    
elseif Var.Speller == 3
    
    for i = 1:sum(cell2mat(Var.NumTrials))
        Var.SC{i} = repmat(7:12,1,Var.NR);
        Var.SR{i} = repmat(1:6,Var.NC,1);
        Var.SR{i} = Var.SR{i}(:)';
    end
    
elseif Var.Speller == 4
    %Changed from this 6/8/17
    %     for i = 1:sum(cell2mat(NumTrials))
    %         SC{i} = repmat(7:12,1,NR);
    %         SR{i} = repmat(1:6,NC,1);
    %         SR{i} = SR{i}(:)';
    %     end
    
    for i = 1:sum(cell2mat(Var.NumTrials))
        Var.SC{i} = 1:Var.NumStim;
        Var.SR{i} = 1:Var.NumStim;
    end
    Var.GoalLetterF = Var.AudioStim(1,TTSpell(:));
    Var.GoalTrialF = 1:sum(cell2mat(Var.NumTrials));
    Var.ChosenLetterF = State.SelectedStimulus(State.SelectedStimulus>0&circshift(State.SelectedStimulus,[1 0])==0);
end


if Var.Speller == 1
    %Check that the markers saved in the log file are the same as those in the
    %data file
    tStimC = double(nonzeros(State.StimulusCodeUNArt)); tStimC = tStimC(1:Var.StateDuration:end);
    clear StimCode
    k=1;
    scind = 0;
    for r = 1:length(Runind)
        for i = 1:Var.NumTrials{r}
            clear sctemp
            for j = 1:Var.NumSequences{r}
                %             sctemp{j} = num2cell(tStimC(NumStimCodes*NumSequences{r}*(i-1)+NumStimCodes*(j-1)+1:...
                %                 NumStimCodes*NumSequences{r}*(i-1)+NumStimCodes*(j)));
                sctemp{j} = num2cell(tStimC([1:Var.NumStimCodes]+scind));
                scind = scind+Var.NumStimCodes;
            end
            StimCode{k} = cell2mat(cat(2,sctemp{:}));
            k=k+1;
        end
    end
    
    %StimCode should contain the same values as mSeq
    k = 1;
    clear isgood
    for r = 1:length(Runind)
        for i = 1:Var.NumTrials{r}
            isgood(k) = isequal(StimCode{k},mSeq{k});
            k=k+1;
        end
    end
    if sum(isgood) ~= sum(cell2mat(Var.NumTrials))
        disp('SOMETHING IS WRONG')
    end
end


%If removing data(probably wont do) need to also remove the same portion of
%StimCode

% NumSequences = NumSequences-SequencestoRemove;
% StimCode = StimCode(1:NumSequences,:);

disp(['Done with myfile ' num2str(toc)]);



end