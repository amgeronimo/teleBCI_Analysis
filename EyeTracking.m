function [gazeacc, gazevar, gazeinvalid, targetbox] = ...
    EyeTracking(GazeX,GazeY, Data, SC,SR, StimulusCode, ...
    SequencePhase, StimulusType, StateDuration, NumSequences, Speller, fs, c, Figures_On, Adopt_Monitor_Specs)


%%
if Adopt_Monitor_Specs == 1
    MonPos =  [1, 1, 1920 1080;1921, 0, 1920 1080];
else
MonPos = get(groot,'MonitorPositions'); %Get monitor positions
end
AppPos = [c.WindowLeft.NumericValue c.WindowTop.NumericValue ...
    c.WindowWidth.NumericValue c.WindowHeight.NumericValue];
COffset = 1000;  %Offset specified in EyetrackerLogger.cpp in ::OnGazeDataEvent
ymondist = MonPos(2,4)-MonPos(1,4)+MonPos(2,2)-1; %find distance from top of monitor 1 to top of monitor 2
MonPos2 = MonPos(2,:); %Isolate second monitor
AppPos(2) = AppPos(2)+COffset;
MonPos2(2) = COffset-ymondist; 

if (isfield(c,'LogEyeTribe') || isfield(c,'LogEyeTribe_v2'))
    MonPos2(1) = 1; MonPos2(2) = MonPos2(2)-COffset; 
    AppPos(1) = AppPos(1)-MonPos(2,1)+1;  AppPos(2) = AppPos(2)-COffset;
end


% figure
% subplot(221)
% scatter(1:length(GazeX),GazeX,[],parula(length(GazeX)))
% line([1 length(GazeX)],MonPos2(1)*[1 1],'Color','b')
% text(length(GazeX)/2,MonPos2(1)*.9,'Left','HorizontalAlignment','center')
% line([1 length(GazeX)],[MonPos2(1)+MonPos2(3)]*[1 1],'Color','r')
% text(length(GazeX)/2,MonPos2(1)+MonPos2(3)*1.1,'Right','HorizontalAlignment','center')
% ylabel('X-postion');
% subplot(223)
% scatter(1:length(GazeY),GazeY,[],parula(length(GazeY)))
% line([1 length(GazeY)],MonPos2(2)*[1 1],'Color','g')
% text(length(GazeY)/2,MonPos2(2)*1.1,'Top','HorizontalAlignment','center')
% line([1 length(GazeY)],(MonPos2(2)+MonPos2(4))*[1 1],'Color','k')
% text(length(GazeY)/2,(MonPos2(2)+MonPos2(4))*1.1,'Bottom','HorizontalAlignment','center')
% ylabel('Y-position');
% subplot(2,2,2)
% scatter(GazeX,GazeY,[],parula(length(GazeX)))
% line(MonPos2(1)*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','b')
% line((MonPos2(1)+MonPos2(3))*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','r')
% line([MonPos2(1) MonPos2(1)+MonPos2(3)],MonPos2(2)*[1 1],'Color','g')
% line([MonPos2(1) MonPos2(1)+MonPos2(3)],(MonPos2(2)+MonPos2(4))*[1 1],'Color','k')
% axis equal;
% xlabel('X-postion'); ylabel('Y-position')




%% Define Gaze quality for each letter

%Define Beginning and End of trials
tmp = find(SequencePhase == 2);
tmp2 = diff(tmp);
tmp3 = find(tmp2>1);
trialboundaries = [[tmp(1); tmp(tmp3+1)] [tmp(tmp3); tmp(end)]];

%Define targets for each trial
targetstim = double(StimulusCode(logical(StimulusType)));
targetstim = targetstim(1:StateDuration:length(targetstim));
targetstim = [targetstim(1:(NumSequences)*2:length(targetstim))...
    targetstim(2:(NumSequences)*2:length(targetstim))];

TargH = c.TargetHeight.NumericValue;
TargW = c.TargetWidth.NumericValue;
BarH = c.StatusBarSize.NumericValue;

%Find locations of these targets
if Speller == 1 %Checkerboard speller
    %Map the target codes to the stimulus
    clear lettertarg
    for ii = 1:length(SC)
    map = [cellfun(@(x) str2num(x),SR{ii}) cellfun(@(x) str2num(x),SC{ii})];
    lettertarg(ii) = find(sum(targetstim(ii,2)==map | targetstim(ii,1)==map,2)>1);
    end
    
    %Didnt want to figure out how to automatically generate the locations,
    %so did it manually
    letterlocs = [1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8 1 2 3 4 5 6 7 8;
                  1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 5 5 5 5 5 5 5 5];
    
    targlocs = letterlocs(:,lettertarg)';

    
    
    
    temp_dist_w = ((100-(TargW*8))/2)/100*AppPos(3);
    temp_dist_h = (((100-BarH)-(TargH*5))/2)/100*AppPos(4)+BarH/100*AppPos(4);
    targleft = temp_dist_w+AppPos(1);
    targtop = temp_dist_h+AppPos(2);
    pxtarg_h = TargH/100*AppPos(4);
    pxtarg_w = TargW/100*AppPos(3);

elseif Speller == 2 %Single Flash (RSVP) speller.
    lettertarg = min(targetstim,[],2);
    targlocs = ones(length(lettertarg),2);
    
        %These are the dimensions in the case of the RSVP speller, using
        %6.25% width and 10% height
%     pxpercm = 35.5;

temp_dist_w = ((100-TargW)/2)/100*AppPos(3);

temp_dist_h = (((100-BarH)-TargH)/2)/100*AppPos(4)+BarH/100*AppPos(4);
targleft = temp_dist_w+AppPos(1);
targtop = temp_dist_h+AppPos(2);
pxtarg_h = TargH/100*AppPos(4);
pxtarg_w = TargW/100*AppPos(3);
    
elseif Speller == 3 %Grid Speller (default P3 task)
    temp = sort(targetstim,2);
    temp(:,2) = temp(:,2)-6;
    lettertarg = 6*(temp(:,1)-1)+temp(:,2);
    targlocs = fliplr(temp);
    
%       pxpercm = 35.5;
%     targleft = 14.8*pxpercm;
%     targtop = 7.05*pxpercm;
%     cmpertarg = 2.95; %The width(and height) of targets.
 temp_dist = ((100-TargSize)/2)/100*AppPos(4);
    targleft = temp_dist+AppPos(1);
    targtop = temp_dist+AppPos(2);
    pxtarg = TargSize/100*AppPos(4);
end


%%Plot

TargSymb = c.TargetDefinitions.Value(:,2)';

%Invert GazeY data and the Y Monitor positions (also invert target Y
%positions in the loop below)
GazeY = -double(GazeY);
MonPos2([2 4]) = -MonPos2([2 4]);
AppPos([2 4]) = -AppPos([2 4]);

if Figures_On == 1
figure('Name','2D position of targets and gaze')
scatter(GazeX,GazeY,[],[.8 .8 .8],'.'); hold on;
line(MonPos2(1)*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','k')
line((MonPos2(1)+MonPos2(3))*[1 1],[MonPos2(2) MonPos2(2)+MonPos2(4)],'Color','k')
line([MonPos2(1) MonPos2(1)+MonPos2(3)],MonPos2(2)*[1 1],'Color','k')
line([MonPos2(1) MonPos2(1)+MonPos2(3)],(MonPos2(2)+MonPos2(4))*[1 1],'Color','k')
line(AppPos(1)*[1 1],[AppPos(2) AppPos(2)+AppPos(4)],'Color',[.5 .5 .5])
line((AppPos(1)+AppPos(3))*[1 1],[AppPos(2) AppPos(2)+AppPos(4)],'Color',[.5 .5 .5])
line([AppPos(1) AppPos(1)+AppPos(3)],AppPos(2)*[1 1],'Color',[.5 .5 .5])
line([AppPos(1) AppPos(1)+AppPos(3)],(AppPos(2)+AppPos(4))*[1 1],'Color',[.5 .5 .5])
line([AppPos(1) AppPos(1)+AppPos(3)],(AppPos(2)+BarH/100*AppPos(4))*[1 1],'Color',[.5 .5 .5])
axis equal;
xlabel('X-postion'); ylabel('Y-position')
clrs = parula(length(lettertarg));

for j = 1:length(lettertarg)
    templeft = targleft + (targlocs(j,1)-1)*pxtarg_w;
    temptop = targtop + (targlocs(j,2)-1)*pxtarg_h;
    targetbox = [templeft -temptop templeft+pxtarg_w...
         -temptop-pxtarg_h];
     line(targetbox(1)*[1 1],[targetbox(2) targetbox(4)],'Color',clrs(j,:));
line(targetbox(3)*[1 1],[targetbox(2) targetbox(4)],'Color',clrs(j,:));
line([targetbox(1) targetbox(3)],targetbox(2)*[1 1],'Color',clrs(j,:));
line([targetbox(1) targetbox(3)],targetbox(4)*[1 1],'Color',clrs(j,:));
disp(['GAZE UPDATES PER TRIAL -- ' num2str(length(unique([find(diff(GazeX(trialboundaries(j,1):trialboundaries(j,2)))>0);...
    find(diff( GazeY(trialboundaries(j,1):trialboundaries(j,2)))>0)])))])
scatter(GazeX(trialboundaries(j,1):trialboundaries(j,2)),...
    GazeY(trialboundaries(j,1):trialboundaries(j,2)),'.',...
    'MarkerEdgeColor',clrs(j,:));
text(targetbox(1),targetbox(2),TargSymb(lettertarg(j)),'Color',...
    clrs(j,:),'FontSize',18,'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
end
end


%% Plot eye movements vs EOG channels
gX = double(GazeX);

if Figures_On == 1
figure('Name','Eye position with EOG channels'); hold on;
plot(Data');
% plot(diff(Data));
plot(100*(gX-mean(gX))/range(gX),'k','LineWidth',2);
end

%% Plot the periods of eye invalidity

DgX = gX-circshift(gX,1);
%find stretches of diff(gX) that are at least .5 seconds long
DgX2 = DgX'~=0;
DgX2 = [1 DgX2 1];
frontSTR = strfind(DgX2,[1 0]);
backSTR = strfind(DgX2,[0 1]);

repeat0length = backSTR-frontSTR;
abovelimit = find(repeat0length>(fs/4));

invaliddata = zeros(length(DgX2),1);
for i = abovelimit
    %the -1 is to remove the 1 added earlier in DgX2
    invaliddata(frontSTR(i):frontSTR(i)+repeat0length(i)-1)=1;
end

if Figures_On == 1
figure('Name','Gaze X position, its derivative, and time of eye gaze invalidity'); hold on;
plot(gX);
plot(DgX);
plot(invaliddata*5000);
end

%Calculate the percentage of the trial for which the gaze was valid
gazeinvalid = zeros(length(lettertarg),1);
for j = 1:length(lettertarg)
gazeinvalid(j) =  sum(invaliddata(trialboundaries(j,1):trialboundaries(j,2)))...
    /diff(trialboundaries(j,:))*100;
end

%% Gaze Accuracy and Variance
gazeacc = zeros(length(lettertarg),2);
gazevar = zeros(length(lettertarg),2);
targetbox = zeros(length(lettertarg),4);
for j = 1:length(lettertarg)
    templeft = targleft + (targlocs(j,1)-1)*pxtarg_w;
    temptop = targtop + (targlocs(j,2)-1)*pxtarg_h;
    targetbox(j,:) = [templeft -temptop templeft+pxtarg_w...
         -temptop-pxtarg_h];
    gxdata = double(GazeX(trialboundaries(j,1):trialboundaries(j,2)));
    gxdata(logical(invaliddata(trialboundaries(j,1):trialboundaries(j,2)))) = NaN;
    gydata = double(GazeY(trialboundaries(j,1):trialboundaries(j,2)));
    gydata(logical(invaliddata(trialboundaries(j,1):trialboundaries(j,2)))) = NaN;
    gazeacc(j,:) =  [nanmean(gxdata)-(templeft+.5*pxtarg_w);
        nanmean(gydata)-(-temptop-.5*pxtarg_h)];
    gazevar(j,:) = [nanvar(gxdata); nanvar(gydata)];
end
%%
if Figures_On == 1
figure('Name','Gaze Statistics','Position',[100 100 900 400])
subplot(131); bar(gazeacc); title('Mean gaze offset');
ax = gca;
ax.XTickLabel = TargSymb(lettertarg);
subplot(132); bar(gazevar); title('Gaze variance');
ax = gca;
ax.XTickLabel = TargSymb(lettertarg);
subplot(133); bar(gazeinvalid); title('Percent invalid gaze');
ax = gca;
ax.XTickLabel = TargSymb(lettertarg);
end


end
