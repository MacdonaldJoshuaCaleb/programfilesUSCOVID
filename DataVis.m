set(0,'DefaultFigureVisible','on')
close all
cla
StatePops = [4903185,731545,7278717,3017804, 39512223, 5758736, 3565287, 973764, 705749,...
    21477737, 10617423,165768, 1415872, 1787065, 12671821, 6732219, 3155070, 2913314, 4467673,...
    4648794, 1344212, 6045680, 6892503, 9986857, 5639632, 2976149, 6137428, 1068778, 1934408,...
    3080156, 1359711, 8882190, 2096829, 19453561, 10488084, 762062,56882, 11689100, 3956971,...
    4217737,12801989,3193694,1059361, 5148714, 884659, 6829174, 28995881, 3205958, 623989,106977,...
    8535519, 7614893,1792147,5822434, 578759]';

load('Data3.mat')
load('USparams.mat')
load('Dfits.mat')
Deaths = Dfits(end,:)';
Data = Data3;
R0s = Data(:,1);
psis = Data(:,2);
alphas = Data(:,3);
captures = Data(:,5);
infecteds = Data(:,7);
F0BS = Data(:,11);
Mortality = Data(:,12);
%Sizes = Data(:,13);

scatter(R0s,psis,30,Mortality,'filled') 
set(gca,'FontSize',16)

% draw the scatter plot
ax = gca;
%ax.XDir = 'reverse';
%view(-31,14)
ylabel('\psi')
xlabel('R_0')
xlim([min(R0s),max(R0s)])
%colormap(flipud(stoplight(64)))
colormap(jet(16))
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Mortality (% population)';


figure
scatter(alphas,psis,30,Mortality,'filled') 


% draw the scatter plot
ax = gca;
%ax.XDir = 'reverse';
%view(-31,14)
ylabel('\psi')
xlabel('\alpha')
xlim([min(alphas),max(alphas)])
colormap(jet(16))
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Mortality (% population)';
set(gca,'FontSize',16)
% map = stoplight(length)
%
%       Returns a custom green-yellow-red colormap by taking only the first
%       third of what would otherwise be MATLAB's standard hsv colormap,
%       and making it slightly darker.
%
% In:
%       length - the length of the colormap.
%
% Out:
%       map - the generated colormap.
% 
%                       Copyright 2017 Laurens R Krol
%                       Team PhyPA, Biological Psychology and Neuroergonomics,
%                       Berlin Institute of Technology
% 2017-03-13 First version
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.





figure
hold on
scatter(R0s,psis,30,alphas,'filled')    % draw the scatter plot
hold off

ax = gca;
%ax.XDir = 'reverse';
%view(-31,14)
xlabel('R_0')
ylabel('\psi')
colormap(jet(16))
xlim([min(R0s),max(R0s)])
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Lockdown Fatigue Measure';
set(gca,'FontSize',16)
figure
hold on

scatter(R0s,psis,30,F0BS,'filled')
hold off

% draw the scatter plot
ax = gca;
%ax.XDir = 'reverse';
%view(-31,14)
xlabel('R_0')
ylabel('\psi')
colormap(jet(16))
xlim([min(R0s),max(R0s)])
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'True Cum. Case Estimate (% pop.)';
set(gca,'FontSize',16)
figure
hold on

scatter(alphas,psis,30,F0BS,'filled')
hold off

% draw the scatter plot
ax = gca;
%ax.XDir = 'reverse';
%view(-31,14)
xlabel('\alpha')
ylabel('\psi')
colormap(jet(16))
xlim([min(alphas),max(alphas)])
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'True Cum. Case Estimate (% pop.)';
set(gca,'FontSize',16)
figure
hold on

scatter(alphas,psis,30,R0s,'filled')
hold off

% draw the scatter plot
ax = gca;
%ax.XDir = 'reverse';
%view(-31,14)
xlabel('\alpha')
ylabel('\psi')
colormap(jet(16))
xlim([min(alphas),max(alphas)])
cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'Reproduction Number';
set(gca,'FontSize',16)
load('R0s.mat')
VarRepos = R0s;
for j = 1:length(VarRepos(:,1))
    for k = 1:55
        if VarRepos(j,k) == 0
            VarRepos(j,k) = NaN;
        end
    end
end
    
figure
hold on
h=fill([77,130,130,77],[0,0,6,6],'k','LineStyle','none');
    h.FaceAlpha=0.1;
plot(0:1:130,ones(1,131),'k-.','linewidth',1.25)
for j = 1:55
    plot(0:1:length(VarRepos(:,1))-1,VarRepos(:,j))
end
hold off
xlim([0 130])
xticks(round(linspace(0,length(VarRepos(:,1))-1,10)))
xlabel('Days Since First Reported Case up to May 31')
ylabel('Time Variable Reproduction Number')
legend({'Stay at Home Order(s) in Effect'},'FontSize',10)
set(gca,'FontSize',16)
load('Cfits.mat')
load('TCTfits.mat')
load('rhos.mat')
load('Tests.mat')
for j = 1:55
    for k = 1:length(rhos(:,1))
    Temp(k,j) = rhos(k,j);
    if Temp(k,j) == 0
        Temp(k,j) = NaN;
    end
    if Temp(k,j) >= 1
        Temp(k,j) = NaN;
    end
    end
end

for j = 1:55
    for k = 1:length(Tests(:,1))
    Temp2(k,j) = 100*(Tests(k,j)/StatePops(j));
    if Temp2(k,j) == 0
       Temp2(k,j) = NaN;
    end
%     if Temp2(k,j) >= 100
%         Temp2(k,j) = NaN;
%     end
    end
end

figure
hold on
h=fill([77,130,130,77],[0,0,max(max(rhos))*1.05,max(max(rhos))*1.05],'k','LineStyle','none');
    h.FaceAlpha=0.1;
for j = 1:55
    plot(0:1:length(VarRepos(:,1))-1,Temp(:,j))
end
hold off
xlim([0 130])
ylim([.95*.03,max(max(rhos))*1.05])
yticks(round(linspace(.05,max(max(rhos)),11),2))
xticks(round(linspace(0,length(VarRepos(:,1))-1,10)))
xlabel('Days Since First Reported Case up to May 31')
ylabel('Ascertainment Rate')
legend({'Stay at Home Order(s) in Effect'},'FontSize',10)
set(gca,'FontSize',16)



format long

figure
hold on
h=fill([77,130,130,77],[0,0,.5,.5],'k','LineStyle','none');
    h.FaceAlpha=0.1;
for j = 1:55
    plot(0:1:length(Tests(:,1))-1,Temp2(:,j))
end
hold off
set(gca,'FontSize',16)
legend({'Stay at Home Order(s) in Effect'},'FontSize',10)
xticks(round(linspace(0,length(VarRepos(:,1))-1,10)))
ylabel('Daily Tests (% State Population)')
xlabel('Days Since First Reproted Case up to May 31')
xlim([0 130])

figure
hold on
%h=fill([77,130,130,77],[0,0,.5,.5],'k','LineStyle','none');
%    h.FaceAlpha=0.1;
for j = 1:55
    plot(Temp2(:,j),Temp(:,j))
end
hold off
set(gca,'FontSize',16)
%legend({'Stay at Home Order(s) in Effect'},'FontSize',10)
%xticks(round(linspace(0,length(VarRepos(:,1))-1,10)))
ylabel('Ascertainment Rate')
xlabel('Daily tests (% State Population)')
