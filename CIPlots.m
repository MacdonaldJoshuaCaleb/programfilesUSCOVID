close all
figure
hold on
%plot(0:1:130,rhoLower,'b',0:1:130,rhoUpper,'b')
plot(0:1:130,rho,'b','linewidth',2)
x2 = [0:1:130,fliplr(0:1:130)];
inBetween = [rhoLower,fliplr(rhoUpper)];
h = fill(x2, inBetween,'b')
set(h,'facealpha',.1)
xlim([0 130])
ylim([0 max(rhoUpper)])
xlabel('Days Since First Reported Case up to May 31')
ylabel('Capture Rate')
legend({'fit \rho(t)','95% Confidence Interval'},'Location','NorthWest')
set(gca,'FontSize',16)
hold off


figure
hold on
%plot(0:1:130,rhoLower,'b',0:1:130,rhoUpper,'b')
plot(0:1:130,R0(1:131),'b','linewidth',2)
x2 = [0:1:130,fliplr(0:1:130)];
inBetween = [R0Lower,fliplr(R0Upper)];
h = fill(x2, inBetween,'b')
set(h,'facealpha',.1)
plot(0:1:130,ones(131,1),'k-.')
xlim([0 130])
ylim([.95*min(R0Lower) 1.05*max(R0Upper)])
xlabel('Days Since First Reported Case up to May 31')
ylabel('Reproduction Number')
legend({'fit R_e','95% Confidence Interval'},'Location','NorthEast')
set(gca,'FontSize',16)
hold off

figure
hold on
%plot(0:1:130,rhoLower,'b',0:1:130,rhoUpper,'b')
plot(0:1:130,(CasesData(1:131)/N)*100,'r*')
plot(0:1:130,(Cfit/N)*100,'b','linewidth',2)
x2 = [0:1:130,fliplr(0:1:130)];
inBetween = [(CasesLower/N)*100,fliplr((CasesUpper/N)*100)];
h = fill(x2, inBetween,'c')
set(h,'facealpha',.1)
plot(0:1:130,(URCfit/N)*100,'b--','linewidth',2)
%x3 = [0:1:130,fliplr(0:1:130)];
inBetween2 = [(URCLower/N)*100,fliplr((URCUpper/N)*100)];
h2 = fill(x2, inBetween2,'b')
set(h2,'facealpha',.1)
plot(0:1:130,((URCfit+Cfit)/N)*100,'k-.','linewidth',2)
%x3 = [0:1:130,fliplr(0:1:130)];
inBetween3 = [((URCLower+CasesLower)/N)*100,fliplr(((CasesLower+URCUpper)/N)*100)];
h3 = fill(x2, inBetween3,'k')
set(h3,'facealpha',.1)
xlim([0 130])
ylim([0 100*1.05*(max(URCUpper)+max(CasesUpper))/N])
xlabel('Days Since First Reported Case up to May 31')
ylabel('% population')
legend({'Reported Case Data','fit Reported Cases','95% CI Reported Cases','fit Unreported Cases','95% CI Unreported Cases','fit True Cases','95% CI True Cases'},'Location','NorthWest','FontSize',10)
set(gca,'FontSize',16)
hold off

figure
hold on
plot(1:1:130,diff(CasesData(1:131)),'r*')
plot(1:1:130,diff(URCfit+Cfit),'k-.','linewidth',2)
x3 = [1:1:130,fliplr(1:1:130)];
inBetween4 = [diff(((URCLower+CasesLower))),fliplr(diff(((URCUpper+CasesUpper))))];
h4 = fill(x3, inBetween4,'k')
set(h4,'facealpha',.1)
plot(1:1:130,diff(Cfit),'b','linewidth',2)
%x3 = [1:1:130,fliplr(1:1:130)];
inBetween5 = [diff(((CasesLower))),fliplr(diff(((CasesUpper))))];
h5 = fill(x3, inBetween5,'c')
set(h5,'facealpha',.1)
xlim([1 130])
ylim([0 1.2*max(diff(URCUpper + CasesUpper))])
legend({'Inferred Daily Case Data','True Daily Case Fit','95% CI True Daily Cases','Daily Reported Case Fit','95% CI Daily Reported Cases'},'Location','NorthWest','FontSize',10)
xlabel('Days Since First Reported Case up to May 31')
ylabel('Count')
set(gca,'FontSize',16)
hold off

