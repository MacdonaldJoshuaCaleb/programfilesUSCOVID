function []= DataSortCleanCasesDeaths
clear
close all
% stop figures from displaying in matlab
set(0,'DefaultFigureVisible','on')
% import clean and shape data to appropriate shape for automated fitting
StatesData = readtable('us-states.csv');
States = table2array(StatesData(:,2));
States = string(States);
StateNames = States(end-54:end);
RawData1 = table2array(StatesData(:,4)); % case data
RawData2 = table2array(StatesData(:,5)); % death data

% data is recorded for 55 states and territories, first case appeared in 
% Washington Sate 55 days ago, DF is solution bin for successfully shaped
% data
count = 0;
for j = 1:size(States)
    if States(j) == "Washington"
        count = count + 1;
    end
end
DF1 = zeros(count,55); % for cases
DF2 = zeros(count,55); % for deaths

%bin for parameter fittings
for k = 1:55
    for j = 1:size(RawData1)
        if States(j) == StateNames(k)
            CData(j) = RawData1(j);
            DData(j) = RawData2(j);
        end
    end
    
CData2 = CData(CData>0);
DData2 = DData(DData>0);
% put data in appropraite location for each state's given row
C1 = zeros(count,1);
D1 = zeros(count,1);
C1(count+1-length(CData2):end) = CData2;
D1(count+1-length(DData2):end) = DData2;
DF1(:,k) = C1;
DF2(:,k) = D1;
clear CData;
clear DData;
end
TestingData = readtable('TestingData.csv');

States2 = table2array(TestingData(:,2));
States2 = string(States2);
counts = zeros(1,length(StateNames));

for j = 1:length(States2)
    for k = 1:length(StateNames)
        Name = StateNames(k);
        if States2(j) == Name
            counts(k) = counts(k)+1;
        end
    end
end
count = max(counts);
DF3 = zeros(count,55);
DF4 = zeros(count,55);
RawData3 = table2array(TestingData(:,4));
RawData4 = table2array(TestingData(:,3));
for k = 1:55
    TData = -ones(1,length(RawData3));
    TCumData = zeros(1,length(RawData3));
    for j = 1:size(RawData3)
        if States2(j) == StateNames(k)
            TData(j) = RawData3(j);
            TCumData(j) = RawData4(j);
        end
    end
    fprintf('State Name is: %s\n',StateNames(k))

TData2 = TData(TData ~=-1);
TCumData2 = TCumData(TCumData >0);
% put data in appropraite location for each state's given row
T1 = zeros(count,1);
TCum1 = zeros(count,1);
T1(count+1-length(TData2):end) = TData2;
TCum1(count+1-length(TCumData2):end) = TCumData2;
DF3(:,k) = T1;
DF4(:,k) = TCum1;
clear TData;
end
save('CasesByState.mat','DF1')
save('DeathsByState.mat','DF2')
save('TestingByState.mat','DF3')
save('CumTestsByState.mat','DF4')
end