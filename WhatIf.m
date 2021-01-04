function [] = WhatIf
close all
set(0,'DefaultFigureVisible','on')
paramfit = [0.286338928636279,258.887865817671,0.0240070128944282,0.00861721915553062,0.358215252628226,271068.991836735,50.3949836983015,12.0000000000000,829407.113710649,1.75000000000000,3.50535973190290,0.0367698558305923,0.0990264973789266,0.0101613689219439];
points = 60;
psis = linspace(.25*paramfit(2),1.75*paramfit(2),points);
alphas = linspace(.005,1/7,points);
betas = [.75*paramfit(1), paramfit(1), 1.25*paramfit(1)];
[X Y] = meshgrid(psis,alphas);

SHE = days(datetime(2020,05,31)-datetime(2020,01,21));
StatePops = [4903185,731545,7278717,3017804, 39512223, 5758736, 3565287, 973764, 705749,...
    21477737, 10617423,165768, 1415872, 1787065, 12671821, 6732219, 3155070, 2913314, 4467673,...
    4648794, 1344212, 6045680, 6892503, 9986857, 5639632, 2976149, 6137428, 1068778, 1934408,...
    3080156, 1359711, 8882190, 2096829, 19453561, 10488084, 762062,56882, 11689100, 3956971,...
    4217737,12801989,3193694,1059361, 5148714, 884659, 6829174, 28995881, 3205958, 623989,106977,...
    8535519, 7614893,1792147,5822434, 578759];
 T = 7.5;
 N = sum(StatePops);
 k0=.005;
deathdist = @(tv,ag,bg) gamcdf(tv,ag,bg);
infdist = @(tv,Ti) expcdf(tv,1/Ti);
tests = [0.460317460317460,0.489795918367347,0.519274376417234,0.548752834467120,0.578231292517007,0.607709750566893,0.637188208616780,0.666666666666667,0.659863945578231,0.653061224489796,0.646258503401361,0.639455782312925,0.632653061224490,0.625850340136054,0.619047619047619,0.544217687074830,0.469387755102041,0.394557823129252,0.319727891156463,0.244897959183673,0.170068027210884,0.0952380952380952,0.136054421768708,0.176870748299320,0.217687074829932,0.258503401360544,0.299319727891156,0.340136054421769,0.380952380952381,129.895691609977,259.410430839002,388.925170068027,518.439909297052,647.954648526077,777.469387755102,906.984126984127,962.251700680272,1017.51927437642,1072.78684807256,1128.05442176871,1183.32199546485,1238.58956916100,1293.85714285714,2593.53061224490,3893.20408163265,5192.87755102041,6492.55102040816,7792.22448979592,9091.89795918367,10391.5714285714,15495.6938775510,20599.8163265306,25703.9387755102,30808.0612244898,35912.1836734694,41016.3061224490,46120.4285714286,54097.1224489796,62073.8163265306,70050.5102040816,78027.2040816327,86003.8979591837,93980.5918367347,101957.285714286,108420.204081633,114883.122448980,121346.040816327,127808.959183673,134271.877551020,140734.795918367,147197.714285714,147751.183673469,148304.653061225,148858.122448980,149411.591836735,149965.061224490,150518.530612245,151072,151668.755102041,152265.510204082,152862.265306122,153459.020408163,154055.775510204,154652.530612245,155249.285714286,166968.612244898,178687.938775510,190407.265306122,202126.591836735,213845.918367347,225565.244897959,237284.571428571,239936.571428571,242588.571428571,245240.571428571,247892.571428571,250544.571428571,253196.571428571,255848.571428571,263339.918367347,270831.265306123,278322.612244898,285813.959183674,293305.306122449,300796.653061224,308288,317228.387755102,326168.775510204,335109.163265306,344049.551020408,352989.938775510,361930.326530612,370870.714285714,375509.836734694,380148.959183674,384788.081632653,389427.204081633,394066.326530612,398705.448979592,403344.571428571,406232.714285714,409120.857142857,412009,414897.142857143,417785.285714286,420673.428571429,423561.571428571,430616.591836735,437671.612244898,444726.632653061,451781.653061225]; 
    function dudt=theModel(t,u,x)
    % u is compartments
    % x is parameter vector, 
    %x = [betag,psig,alphag,phig,Tg,Tqg,xig,kg,Ag]
    beta = x(1);  
    psi = x(2);    
    alpha = x(3); 
    phi = 0;   
    %T = x(9);     
  
    betaq = 0;
    nu = 0;
    phiq = 0;
   
%     tests=interp1(tD0,avs,t);
%     rho = k*tests./(A+tests)+k0;
    dudt = [
        -(beta/N*u(3))*u(1)*(1+psi) + alpha*u(2);
        psi*(beta/N*u(3))*u(1) - u(2)*(alpha);
        (beta/N*u(3))*u(1) - (1/T)*u(3)
        ];
  end
function z = objFxn2(x,tD)
        s0 = N;
        sq0 = x(9);
        k = x(5);
        A = x(6);
        z1 = zeros(length(tD),1);
        z2 = zeros(length(tD),1);
        z3 = zeros(length(tD),1);
        %tests=movmean(NewTests,3);
        rho = k*((tests/N))./(A/N+((tests/N)))+k0;
        mdd= zeros(length(tD),1);
        mid= zeros(length(tD),1);
        z1(1)=1;
        z2(1)=0;
        z3(1)=x(7)-1;
%         I0 = Cases(1)/rho(1);
        I0 = x(7);
        ag = x(8);
        bg = 21/ag;
        Iq0 = 0;
     
        %T = x(9);     
      
        u0 = [s0-sq0-I0,sq0,I0];
        xi = x(4);
        [t,Sol] = ode45(@(t,u)theModel(t,u,x),tD,u0);
        infs = Sol(:,3);
        infsq =  zeros(length(infs),1);
       % defs = Sol(:,5) + Sol(:,6);
       % z1= Sol(:,5); 
       % mass(1) = (gamcdf(1.5,5.1,.86)-gamcdf(0,5.1,.86))+(gamcdf(1.5,17.8,.45)-gamcdf(0,17.8,.45));
        for j = 2:length(tD)
%             mass(j) =(gamcdf(j+.5,5.1,.86)-gamcdf(j-.5,5.1,.86))+(gamcdf(j+.5,17.8,.45)-gamcdf(j-.5,17.8,.45));
%             defs(j) = mass(j)*defs(j);
            mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1);
            mdd(j-1)=(deathdist((2:j),ag,bg)-deathdist((1:j-1),ag,bg))*(infs(j-1:-1:1)+infsq(j-1:-1:1));
            z2(j) = z2(j-1)+mdd(j-1);
            if j <= SHE
            z1(j) = z1(j-1)+rho(j-1)*mid(j-1);
            z3(j) = z3(j-1)+(1-rho(j-1))*mid(j-1);
            end
            % get rid z2(j-1), z1(j-1)
        end
%         for j = 1:length(tD)
%             z1(j) = sum(infs(1:j));    
%         end
        z2 = z2*xi;
        z = [z1(1:SHE)',z2(1:SHE+21)',z3(1:SHE)'];
end
tD = 0:1:SHE+21;
fits = objFxn2(paramfit,tD);
Cfit = fits(1:SHE);
Dfit = fits(SHE+1:2*SHE+21);
URCfit = fits(2*SHE+22:end);
save('USCfitCurr.mat','Cfit')
save('USDfitCurr.mat','Dfit')
save('USURCfitCurr.mat','URCfit')
fprintf('Stop Here\n')
MaxDailyCases1 = zeros(points,points);
MaxDailyCases2 = zeros(points,points);
MaxDailyCases3 = zeros(points,points);
CumCase1 = zeros(points,points);
CumCase2 = zeros(points,points);
CumCase3 = zeros(points,points);
CumDeaths1 = zeros(points,points);
CumDeaths2 = zeros(points,points);
CumDeaths3 = zeros(points,points);
DateMax1 = zeros(points,points);
DateMax2 = zeros(points,points);
DateMax3 = zeros(points,points);
% for p = 1:3
%     for q = 1:points
%         for r = 1:points
%             if p == 1
%             param = [betas(p),X(q,r),Y(q,r),paramfit(4:9)];
%             fits = objFxn2(param,tD);
%             RCfit = fits(1:SHE);
%  %           Dfit = fits(SHE+1:2*SHE+21);
%             URCfit = fits(2*SHE+22:end);
%             Cfit = RCfit + URCfit;
%             MaxDailyCases1(q,r) = (max(diff(Cfit))/N)*100;
%             [M, I] = max(diff(Cfit));
%             DateMax1(q,r) = I;
%             CumCase1(q,r) = (Cfit(end)/N)*100;
% %             CumDeaths1(q,r) = (Dfit(end)/N)*100;
%             end
%             if p == 2
%             param = [betas(p),X(q,r),Y(q,r),paramfit(4:9)];
%             fits = objFxn2(param,tD);
%             RCfit = fits(1:SHE);
% %             Dfit = fits(SHE+1:2*SHE+21);
%             URCfit = fits(2*SHE+22:end);
%             Cfit = RCfit + URCfit;
%             MaxDailyCases2(q,r) = (max(diff(Cfit))/N)*100;
%             [M, I] = max(diff(Cfit));
%             DateMax2(q,r) = I;
%             CumCase2(q,r) = (Cfit(end)/N)*100;
% %             CumDeaths2(q,r) = (Dfit(end)/N)*100;
%             end
%             if p == 3
%             param = [betas(p),X(q,r),Y(q,r),paramfit(4:9)];
%             fits = objFxn2(param,tD);
%             RCfit = fits(1:SHE);
% %             Dfit = fits(SHE+1:2*SHE+21);
%             URCfit = fits(2*SHE+22:end);
%             Cfit = RCfit + URCfit;
%             MaxDailyCases3(q,r) = (max(diff(Cfit))/N)*100;
%             [M, I] = max(diff(Cfit));
%             DateMax3(q,r) = I;
%             CumCase3(q,r) = (Cfit(end)/N)*100;
% %             CumDeaths3(q,r) = (Dfit(end)/N)*100;
%             end
%         end
%     end
% end
% 
% 
% 
% 
% figure
% tiledlayout(3,1)
% ax1 = nexttile;
% contourf(X,Y,MaxDailyCases3);
% title('R_0=2.69')
% %colormap(ax1,winter)
% colorbar
% ax2 = nexttile; 
% hold on
% contourf(X,Y,MaxDailyCases2);
% plot(paramfit(2),paramfit(3),'r*','linewidth',2)
% hold off
% title('R_0=2.15')
% ylabel('\alpha','FontSize',16)
% %colormap(ax2,winter)
% colorbar
% ax3 = nexttile;
% contourf(X,Y,MaxDailyCases1);
% title('R_0=1.61')
% xlabel('\psi','FontSize',16)
% %colormap(ax3,winter)
% colorbar
% sgtitle('Peak Daily Case Incidence (% population)','FontSize',16)
% 
% figure
% tiledlayout(3,1)
% ax1 = nexttile;
% contourf(X,Y,CumCase3);
% title('R_0=2.69')
% %colormap(ax1,winter)
% colorbar
% ax2 = nexttile; 
% hold on
% contourf(X,Y,CumCase2);
% plot(paramfit(2),paramfit(3),'r*','linewidth',2)
% hold off
% title('R_0=2.15')
% ylabel('\alpha','FontSize',16)
% %colormap(ax2,winter)
% colorbar
% ax3 = nexttile;
% contourf(X,Y,CumCase1);
% title('R_0=1.61')
% xlabel('\psi','FontSize',16)
% %colormap(ax3,winter)
% colorbar
% sgtitle('Cumlative Case Total (% population)','FontSize',16)
% 
% % figure
% % tiledlayout(3,1)
% % ax1 = nexttile;
% % contourf(X,Y,CumDeaths3);
% % title('R_0=2.69')
% % %colormap(ax1,winter)
% % colorbar
% % ax2 = nexttile; 
% % hold on
% % contourf(X,Y,CumDeaths2);
% % plot(paramfit(2),paramfit(3),'r*','linewidth',2)
% % hold off
% % title('R_0=2.15')
% % ylabel('\alpha','FontSize',16)
% % %colormap(ax2,winter)
% % colorbar
% % ax3 = nexttile;
% % contourf(X,Y,CumDeaths1);
% % title('R_0=1.61')
% % xlabel('\psi','FontSize',16)
% % %colormap(ax3,winter)
% % colorbar
% % sgtitle('Cumlative Death Total (% population)','FontSize',16)
% 
% figure
% tiledlayout(3,1)
% ax1 = nexttile;
% contourf(X,Y,DateMax3);
% title('R_0=2.69')
% 
% colorbar
% ax2 = nexttile; 
% hold on
% contourf(X,Y,DateMax2);
% plot(paramfit(2),paramfit(3),'r*','linewidth',2)
% hold off
% title('R_0=2.15')
% ylabel('\alpha','FontSize',16)
% 
% colorbar
% ax3 = nexttile;
% contourf(X,Y,DateMax1);
% title('R_0=1.61')
% xlabel('\psi','FontSize',16)
% 
% colorbar
% sgtitle({'Time to Peak Daily Case Incidence (Days)'},'FontSize',16)

timefits = zeros(length(betas),130);
for kk = 1:length(betas)
param = [betas(kk),paramfit(2),paramfit(3),paramfit(4:9)];
fits = objFxn2(param,tD);
timefits(kk,:) = diff(fits(1:SHE)+fits(2*SHE+22:end));
u0 = [N-param(9)-param(7),param(9),param(7)];
[t,Sol] = ode45(@(t,u)theModel(t,u,param),tD,u0);
quars(kk,:) = Sol(1:SHE,2);
clear Sol
end
M = .2;
fig = figure;
   red = [0 0 0];
   blue = [0 0 0];
   left_color = blue;
   right_color = red;
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(1:1:SHE-1, (timefits(1,:)/N)*100,'k-', 1:1:SHE-1, (timefits(2,:)/N)*100,'r-*', 1:1:SHE-1, (timefits(3,:)/N)*100,'b-','linewidth',2)
xlabel('Days Since First Reported Case up to May 31')
ylabel('Daily Case Incidence (% pop)')

ylim([0, 1.05*M])
xlim([1, SHE-1])
%title(Names(its(ss)))
yyaxis right
plot(1:1:SHE-1, (quars(1,2:SHE)/N)*100,'k--',1:1:SHE-1, (quars(2,2:SHE)/N)*100,'r--',1:1:SHE-1, (quars(3,2:SHE)/N)*100,'b--','linewidth',1)
ylabel('Quarantined  Individuals (% pop)')
str1 = strcat('R_0 =',' ',num2str(round(betas(1)*7.5,2)));
str2 = strcat('R_0 =',' ',num2str(round(betas(2)*7.5,2)),' (fit value)');
str3 = strcat('R_0 =',' ',num2str(round(betas(3)*7.5,2)));
legend({str1,str2,str3},'FontSize',10,'Location','NorthWest')  
set(gca,'FontSize',16)

timefits = zeros(length(betas),130);
psisT = [.75*paramfit(2),paramfit(2),1.75*paramfit(2)];
for kk = 1:length(betas)
param = [paramfit(1),psisT(kk),paramfit(3),paramfit(4:9)];
fits = objFxn2(param,tD);
timefits(kk,:) = diff(fits(1:SHE)+fits(2*SHE+22:end));
u0 = [N-param(9)-param(7),param(9),param(7)];
[t,Sol] = ode45(@(t,u)theModel(t,u,param),tD,u0);
quars(kk,:) = Sol(1:SHE,2);
clear Sol
end

fig = figure;
   red = [0 0 0];
   blue = [0 0 0];
   left_color = blue;
   right_color = red;
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(1:1:SHE-1, (timefits(1,:)/N)*100,'b-', 1:1:SHE-1, (timefits(2,:)/N)*100,'r-*', 1:1:SHE-1, (timefits(3,:)/N)*100,'k-','linewidth',2)
xlabel('Days Since First Reported Case up to May 31')
ylabel('Daily Case Incidence (% pop)')

ylim([0, 1.05*M])
xlim([1, SHE-1])
%title(Names(its(ss)))
yyaxis right
plot(1:1:SHE-1, (quars(1,2:SHE)/N)*100,'b--',1:1:SHE-1, (quars(2,2:SHE)/N)*100,'r--',1:1:SHE-1, (quars(3,2:SHE)/N)*100,'k--','linewidth',1)
ylabel('Quarantined  Individuals (% pop)')
str1 = strcat('\psi =',' ',num2str(round(psisT(1),1)));
str2 = strcat('\psi =',' ',num2str(round(psisT(2),1)),' (fit value)');
str3 = strcat('\psi =',' ',num2str(round(psisT(3),1)));
legend({str1,str2,str3},'FontSize',10,'Location','NorthWest')  
set(gca,'FontSize',16)

timefits = zeros(length(betas),130);
alphasT = [.005,paramfit(3),1/7];
for kk = 1:length(betas)
param = [paramfit(1),paramfit(2),alphasT(kk),paramfit(4:9)];
fits = objFxn2(param,tD);
timefits(kk,:) = diff(fits(1:SHE)+fits(2*SHE+22:end));
u0 = [N-param(9)-param(7),param(9),param(7)];
[t,Sol] = ode45(@(t,u)theModel(t,u,param),tD,u0);
quars(kk,:) = Sol(1:SHE,2);
clear Sol
end

fig = figure;
   red = [0 0 0];
   blue = [0 0 0];
   left_color = blue;
   right_color = red;
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(1:1:SHE-1, (timefits(1,:)/N)*100,'k-', 1:1:SHE-1, (timefits(2,:)/N)*100,'r-*', 1:1:SHE-1, (timefits(3,:)/N)*100,'b-','linewidth',2)
xlabel('Days Since First Reported Case up to May 31')
ylabel('Daily Case Incidence (% pop)')

ylim([0, 1.05*M])
xlim([1, SHE-1])
%title(Names(its(ss)))
yyaxis right
plot(1:1:SHE-1, (quars(1,2:SHE)/N)*100,'k--',1:1:SHE-1, (quars(2,2:SHE)/N)*100,'r--',1:1:SHE-1, (quars(3,2:SHE)/N)*100,'b--','linewidth',1)
ylabel('Quarantined  Individuals (% pop)')
str1 = strcat('\alpha =',' ',num2str(round(1/alphasT(1),1)));
str2 = strcat('\alpha =',' ',num2str(round(1/alphasT(2),1)),' (fit value)');
str3 = strcat('\alpha =',' ',num2str(round(1/alphasT(3),1)));
legend({str1,str2,str3},'FontSize',10,'Location','NorthWest')  
title('United States')
set(gca,'FontSize',16)
end
