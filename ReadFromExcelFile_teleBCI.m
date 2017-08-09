function[PItems, PData] = ReadFromExcelFile_teleBCI(filename)

hExcel = actxserver('Excel.Application');
hExcel.visible = 1; % If you want Excel visible.
hExcel.DisplayAlerts = false; % Avoid excel warning popups


cdp = cd;
% [~, sl] = regexp(cdp,'\Users\');
% se = regexp(cdp,'Box Sync');
% usrnm = cdp(sl+1:se-2);
% pathToFile = ['C:\Users\' usrnm '\Box Sync\Hershey_2016\ALSA teleBCI\Analysis'];
%

Wkbk = hExcel.Workbooks.Open(fullfile(cdp,[filename '.csv'])); % Opens Excel file
Patient = Wkbk.Sheets.Item(1); % Get the Patient sheet
PRange = Patient.UsedRange;

PData = PRange.Value(2:end,:);

PItems = PRange.Value(1,:);


%     for pp = 1:length(PItems)
%         [row col] = find(strcmp(PRange.Value(1,:),PItems{pp})==1);
%         PData{pp}= PatDat(:,col);
%     end


%    %Sheets.Item(1).Delete; % Deletes the first sheet
%Wkbk.ActiveSheet.Cells.Item(1).EntireRow.Insert; Assuming you want to insert a row and not a column
%    Wkbk.ActiveSheet.Cells.Range('A1').Value = filenames{i}; % Insert Filename
%    Wkbk.Save; Save xls file
%    Wkbk.SaveAs('C:\MatlabWork\Answers\Test.csv', 6); Save as .csv (file format for csv enum is 6)
Wkbk.Close;


hExcel.Quit;
hExcel.delete;
end