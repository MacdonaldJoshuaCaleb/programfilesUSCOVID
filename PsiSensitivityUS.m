StatePops = [4903185,731545,7278717,3017804, 39512223, 5758736, 3565287, 973764, 705749,...
    21477737, 10617423,165768, 1415872, 1787065, 12671821, 6732219, 3155070, 2913314, 4467673,...
    4648794, 1344212, 6045680, 6892503, 9986857, 5639632, 2976149, 6137428, 1068778, 1934408,...
    3080156, 1359711, 8882190, 2096829, 19453561, 10488084, 762062,56882, 11689100, 3956971,...
    4217737,12801989,3193694,1059361, 5148714, 884659, 6829174, 28995881, 3205958, 623989,106977,...
    8535519, 7614893,1792147,5822434, 578759]';
set(0,'DefaultFigureVisible','on')
close all
load('CasesWhaatIf.mat')
load('CasesWhatIfAlpha.mat')
load('PsisWhatIf.mat')
load('peakCases.mat')
load('peakCasesAlpha.mat')
Table = readtable('PsiSensativity.csv');
Names = table2array(Table(:,1));
Names = string(Names);
CasesWhatifAlpha(:,9) = [];
CasesWhatifAlpha(:,12) = [];
CasesWhatifAlpha(:,37) = [];
CasesWhatifAlpha(:,42) = [];
CasesWhatifAlpha(:,50) = [];

CasesWhatif(:,9) = [];
CasesWhatif(:,12) = [];
CasesWhatif(:,37) = [];
CasesWhatif(:,42) = [];
CasesWhatif(:,50) = [];

PsisWhatIf(:,9) = [];
PsisWhatIf(:,12) = [];
PsisWhatIf(:,37) = [];
PsisWhatIf(:,42) = [];
PsisWhatIf(:,50) = [];

peakCasesAlpha(:,9) = [];
peakCasesAlpha(:,12) = [];
peakCasesAlpha(:,37) = [];
peakCasesAlpha(:,42) = [];
peakCasesAlpha(:,50) = [];

peakCases(:,9) = [];
peakCases(:,12) = [];
peakCases(:,37) = [];
peakCases(:,42) = [];
peakCases(:,50) = [];

SortByR0 = sortrows(Table,2);
R0s = table2array(SortByR0(:,2))*7.5;
SortByCumCases = sortrows(Table,5);
SortByPeakCases = sortrows(Table,6);
IndexR0 = table2array(SortByR0(:,7));
IndexCC = table2array(SortByCumCases(:,7));
IndexPC = table2array(SortByPeakCases(:,7));

figure
tiledlayout(2,1)
ax1 = nexttile;
hold on
cols = jet(50);
colormap(cols);
  cb = colorbar; 
  caxis([min(R0s), max(R0s)]);
 cb.Label.String = 'R_0';

for j = 1:50
    plot(PsisWhatIf(:,IndexR0(j)),CasesWhatif(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
end 
for j = 1:50
    p = plot(PsisWhatIf(25,IndexR0(j)),CasesWhatif(25,IndexR0(j)),'linewidth',2,'Color',cols(j,:));
    p.Marker = '.';
    %p.MarkerEdgeColor = 'k';
    %p.MarkerFaceColor = cols(j,:);
    p.MarkerSize = 40;
end

%for j = 1:50
%    plot(PsisWhatIf(:,IndexR0(j)),CasesWhatif(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
%end


hold off
title('Sensitivity of Cum. Cases to Change in \psi, \alpha, R_0')
xlabel('\psi (log scale)')
ylabel('Cum. Cases (% population)')
xlim([1 6000])
ylim([0 100])
lgd=legend();
lgd.NumColumns = 2;
lgd.FontSize = 12;
lgd.Location = 'NorthWestOutside';
lgd.Title.String = 'States Ordered by Increasing Reproduction Number';
set(gca,'FontSize',16)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax2 = nexttile;
hold on
cols = jet(50);
colormap(cols);
  cb = colorbar; 
  caxis([min(R0s), max(R0s)]);
 cb.Label.String = 'R_0';
 


for j = 1:50
    plot(PsisWhatIf(:,IndexR0(j)),CasesWhatifAlpha(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
end

for j = 1:50
    p = plot(PsisWhatIf(25,IndexR0(j)),CasesWhatif(25,IndexR0(j)),'linewidth',2,'Color',cols(j,:));
    p.Marker = '.';
    %p.MarkerEdgeColor = 'k';
    %p.MarkerFaceColor = cols(j,:);
    p.MarkerSize = 40;
end
hold off
title('Sensitivity of Cum. Cases to Change in \psi, R_0 (\alpha = 0)')
xlabel('\psi (log scale)')
ylabel('Cum. Cases (% population)')
xlim([1 6000])
ylim([0 100])
%lgd=legend();
%lgd.NumColumns = 2;
%lgd.FontSize = 12;
%lgd.Location = 'NorthWestOutside';
%lgd.Title.String = 'States Ordered by Increasing Reproduction Number';
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
tiledlayout(2,1)
ax1 = nexttile;
hold on
cols = jet(50);
colormap(cols);
  cb = colorbar; 
  caxis([min(R0s), max(R0s)]);
 cb.Label.String = 'R_0';

for j = 1:50
    plot(PsisWhatIf(:,IndexR0(j)),peakCases(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
end 
for j = 1:50
    p = plot(PsisWhatIf(25,IndexR0(j)),peakCases(25,IndexR0(j)),'linewidth',2,'Color',cols(j,:));
    p.Marker = '.';
    %p.MarkerEdgeColor = 'k';
    %p.MarkerFaceColor = cols(j,:);
    p.MarkerSize = 40;
end

%for j = 1:50
%    plot(PsisWhatIf(:,IndexR0(j)),CasesWhatif(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
%end


hold off
title('Sensitivity of Peak Cases to Change in \psi, \alpha, R_0')
xlabel('\psi (log scale)')
ylabel('Cum. Cases (% population)')
xlim([1 6000])
ylim([0 max(max(peakCasesAlpha))])
lgd=legend();
lgd.NumColumns = 2;
lgd.FontSize = 12;
lgd.Location = 'NorthWestOutside';
lgd.Title.String = 'States Ordered by Increasing Reproduction Number';
set(gca,'FontSize',16)





ax2 = nexttile;
hold on
cols = jet(50);
colormap(cols);
  cb = colorbar; 
  caxis([min(R0s), max(R0s)]);
 cb.Label.String = 'R_0';
 


for j = 1:50
    plot(PsisWhatIf(:,IndexR0(j)),peakCasesAlpha(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
end

for j = 1:50
    p = plot(PsisWhatIf(25,IndexR0(j)),peakCases(25,IndexR0(j)),'linewidth',2,'Color',cols(j,:));
    p.Marker = '.';
    %p.LineWidth = 2;
    %p.MarkerEdgeColor = 'k';
    %p.MarkerFaceColor = cols(j,:);
    p.MarkerSize = 40;
end

hold off
title('Sensitivity of Peak Cases to Change in \psi, R_0 (\alpha = 0)')
xlabel('\psi (log scale)')
ylabel('Cum. Cases (% population)')
xlim([1 6000])
ylim([0 max(max(peakCasesAlpha))])
%lgd=legend();
%lgd.NumColumns = 2;
%lgd.FontSize = 12;
%lgd.Location = 'NorthWestOutside';
%lgd.Title.String = 'States Ordered by Increasing Reproduction Number';
set(gca,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
tiledlayout(2,1)
ax1 = nexttile;

hold on
cols = jet(50);
colormap(cols);
  cb = colorbar; 
  caxis([min(R0s), max(R0s)]);
 cb.Label.String = 'R_0';
 


for j = 1:50
    plot(PsisWhatIf(:,IndexR0(j)),CasesWhatifAlpha(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
end

for j = 1:50
    p = plot(PsisWhatIf(25,IndexR0(j)),CasesWhatif(25,IndexR0(j)),'linewidth',2,'Color',cols(j,:));
    p.Marker = '.';
    %p.MarkerEdgeColor = 'k';
    %p.MarkerFaceColor = cols(j,:);
    p.MarkerSize = 40;
end
hold off
title('Sensitivity of Cum. Cases to Change in \psi, R_0 (\alpha = 0)')
xlabel('\psi (log scale)')
ylabel('Cum. Cases (% population)')
xlim([1 6000])
ylim([0 100])
%lgd=legend();
%lgd.NumColumns = 2;
%lgd.FontSize = 12;
%lgd.Location = 'NorthWestOutside';
%lgd.Title.String = 'States Ordered by Increasing Reproduction Number';
set(gca,'FontSize',16)
xlim([1 6000])
ylim([0 max(max(CasesWhatifAlpha))])
% lgd=legend();
% lgd.NumColumns = 5;
% lgd.FontSize = 12;
% lgd.Location = 'SouthWestOutside';
% lgd.Title.String = 'States Ordered by Increasing Reproduction Number';
% set(gca,'FontSize',16)





ax2 = nexttile;
hold on
cols = jet(50);
colormap(cols);
  cb = colorbar; 
  caxis([min(R0s), max(R0s)]);
 cb.Label.String = 'R_0';
 


for j = 1:50
    plot(PsisWhatIf(:,IndexR0(j)),peakCasesAlpha(:,IndexR0(j)),'Color',cols(j,:),'linewidth',2,'DisplayName',Names(IndexR0(j)))
end

for j = 1:50
    p = plot(PsisWhatIf(25,IndexR0(j)),peakCases(25,IndexR0(j)),'linewidth',2,'Color',cols(j,:));
    p.Marker = '.';
    %p.LineWidth = 2;
    %p.MarkerEdgeColor = 'k';
    %p.MarkerFaceColor = cols(j,:);
    p.MarkerSize = 40;
end

hold off
title('Sensitivity of Peak Cases to Change in \psi, R_0 (\alpha = 0)')
xlabel('\psi (log scale)')
ylabel('Cum. Cases (% population)')
xlim([1 6000])
ylim([0 max(max(peakCasesAlpha))])
lgd=legend();
lgd.NumColumns = 5;
lgd.FontSize = 12;
lgd.Location = 'NorthWestOutside';
lgd.Title.String = 'States Ordered by Increasing Reproduction Number';
set(gca,'FontSize',16)
set(gca,'FontSize',16)

