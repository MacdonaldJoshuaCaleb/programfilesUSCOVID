function [] = USStateFitsSimplest
warning off;
close all
set(0,'DefaultFigureVisible','off')
% import data
paramfits = readtable('params.csv');
paramfits = table2array(paramfits);
CasesByState = readtable('CasesByState.csv');
Dates = table2array(CasesByState(:,1));
CasesByState = table2array(CasesByState(:,2:end));
DeathsByState = readtable('DeathsByState.csv');
DeathsByState = table2array(DeathsByState(:,2:end));
TestingByState = readtable('CumTestsByState.csv');
TestingByState = table2array(TestingByState(:,2:end));
DatesByState = readtable('ImportantDates.csv');
DatesByState = table2array(DatesByState(:,2:end));
AtHomeByState = readtable('AtHomeByState.csv');
DatesAtHome = table2array(AtHomeByState(:,1));
AtHomeByState = table2array(AtHomeByState(:,2:end));
StateNames=["Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","District of Columbia","Florida","Georgia","Guam","Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Northern Mariana Islands","Ohio","Oklahoma","Oregon","Pennsylvania","Puerto Rico","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virgin Islands","Virginia","Washington","West Virginia","Wisconsin","Wyoming"];
Names = StateNames;
StatePops = [4903185,731545,7278717,3017804, 39512223, 5758736, 3565287, 973764, 705749,...
    21477737, 10617423,165768, 1415872, 1787065, 12671821, 6732219, 3155070, 2913314, 4467673,...
    4648794, 1344212, 6045680, 6892503, 9986857, 5639632, 2976149, 6137428, 1068778, 1934408,...
    3080156, 1359711, 8882190, 2096829, 19453561, 10488084, 762062,56882, 11689100, 3956971,...
    4217737,12801989,3193694,1059361, 5148714, 884659, 6829174, 28995881, 3205958, 623989,106977,...
    8535519, 7614893,1792147,5822434, 578759];

 T = 7.5;
% lsqwt = .05;
% remove Louisiana as already handled 

% stay home start and end
SHSs = DatesByState(1,:);
SHEs = DatesByState(2,:);

% get date of first reported case
for w = 1:length(Names)
    for j = 1:length(Dates)
        if CasesByState(j,w) > 0
            FirstCases(w) = Dates(j);
            break
        end
    end
end

FirstCases = FirstCases';
FirstCases(52)


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

function z = objFxn(x,tD)
        s0 = N;
        sq0 = x(9);
        k = x(5);
        A = x(6);
        z1 = zeros(length(tD),1);
        z2 = zeros(length(tD),1);
        %tests=movmean(NewTests,3);
        tests = interp1(tD0,avs,tD);
        rho = k*((tests/N))./(A/N+((tests/N)))+k0;
        mdd= zeros(length(tD),1);
        mid= zeros(length(tD),1);
        z1(1)=Cases(1);
        z2(1)=Deaths(1);
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
            z1(j) = z1(j-1)+rho(j-1)*mid(j-1);
            % get rid z2(j-1), z1(j-1)
        end
%         for j = 1:length(tD)
%             z1(j) = sum(infs(1:j));    
%         end
        z2 = z2*xi;
        z = [lsqwt*z1(1:SHE)',z2(1:SHE+21)'];
end

rhos = zeros(days(datetime(2020,05,31)-FirstCases(52)),55);
R0s = zeros(days(datetime(2020,05,31)-FirstCases(52)),55);
Cfits = zeros(days(datetime(2020,05,31)-FirstCases(52)),55);
Dfits = zeros(days(datetime(2020,05,31)-FirstCases(52)+21),55);
URCfits = zeros(days(datetime(2020,05,31)-FirstCases(52)),55);
Positives = zeros(days(datetime(2020,05,31)-FirstCases(52)-1),55);
Tests = zeros(days(datetime(2020,05,31)-FirstCases(52)),55);
PsisWhatIf = zeros(25,55);
CasesWhatif = zeros(25,55);
CasesWhatifAlpha = zeros(25,55);
% Fitting Loop
for its = 1:length(Names)
    Cases = CasesByState(:,its);
    Deaths = DeathsByState(:,its);
    NewTests = diff(TestingByState(:,its));
    StayHomeStart = SHSs(its);
    StayHomeEnd = SHEs(its);
    FirstCase = FirstCases(its);
    SHE = days(StayHomeEnd-FirstCase);
    SHS = days(StayHomeStart-FirstCase);
    TF = isnan(SHE);
    if TF == 1
        SHE = 59;
        SHS = 0;
    end
    if SHE <= 59
        SHE = 59;
    end
    SHE = days(datetime(2020,05,31)-FirstCase);
        
    

    
    
    % shape particular state data appropriately 
    Cases = Cases(Cases>0);
    

    
    for w = 1:length(NewTests)
        if NewTests(w) > 0
            break
        end
    end
    NewTests = NewTests(w:end);
    
    
    
    
    if NewTests(1) == 0
    %Cd(1) = 2*Cd(2)-Cd(3);
        NewTests(1) = (NewTests(2))/2;
        NewTests(2)= (NewTests(2))/2;
        if NewTests(2) == 0
            NewTests(1) = NewTests(3)/3;
            NewTests(2) = NewTests(3)/3;
            NewTests(3) = NewTests(3)/3;
        end
    end
    if NewTests(end) == 0
    %Cd(1) = 2*Cd(2)-Cd(3);
        NewTests(end) = (NewTests(end-1))/2;
        NewTests(end-1)= (NewTests(end-1))/2;
        if NewTests(end-1) == 0
            NewTests(end) = NewTests(end-2)/3;
            NewTests(end-1) = NewTests(end-2)/3;
            NewTests(end-2) = NewTests(end-2)/3;
        end
    end
    if NewTests(1) < 0
        NewTests(1) = -NewTests(1);
    end
    ind0=find(NewTests<0);
    if isempty(ind0)==0
        NewTests(ind0)=1/3*NewTests(ind0-1)+1/3*NewTests(ind0+1);
        NewTests(ind0-1)=2/3*NewTests(ind0-1);
        NewTests(ind0+1)=2/3*NewTests(ind0+1);
    end
    clear ind0;
    cc=1;
    while NewTests(cc) < .75*max(NewTests(1:SHE))
        cc = cc+1;
    end
    N = StatePops(its);
    Deaths = Deaths(length(Deaths)-length(Cases)+1:end);
    % get weekly testing averages 
    intv=3;
    Weeks = floor(length(NewTests)/intv);
    avs = zeros(length(Weeks)+1,1);
    for w = 1:Weeks
        avs(w) = mean(NewTests(1+intv*(w-1):intv+intv*(w-1)));
    end
    avs(Weeks+1) = mean(NewTests(Weeks*intv+1:end));
    tW = linspace(0,Weeks,Weeks+1);
    tD0=tW*intv;
    k0=.03;
    lsqwt=max(Deaths(1:SHE+21))/max(Cases(1:SHE));

   
  

    deathdist = @(tv,ag,bg) gamcdf(tv,ag,bg);
    infdist = @(tv,Ti) expcdf(tv,1/Ti);
      
        

    tD = linspace(0,length(Cases)-1,length(Cases));
    Cd = diff(Cases);
    Dd = diff(Deaths);
    tdc = tD(2:end)-1;

    if Cd(1) == 0
        temp = 2;
        while Cd(1) == 0
            for j = 1:temp
                Cd(j) = Cd(temp)/temp;
            end
            temp = temp+1;
        end
    end
    if Cd(end) == 0
        %Cd(1) = 2*Cd(2)-Cd(3);
        Cd(end) = (Cd(end-1))/2;
        Cd(end-1)= (Cd(end-1))/2;
        if Cd(end-1) == 0
            Cd(end) = Cd(end-2)/3;
            Cd(end-1) = Cd(end-2)/3;
            Cd(end-2) = Cd(end-2)/3;
        end
    end
    if Cd(1) < 0
        Cd(1) = -Cd(1);
    end
    Cd(1)
    ind0=find(Cd<0);
    if isempty(ind0)==0
        Cd(ind0)=1/3*Cd(ind0-1)+1/3*Cd(ind0+1);
        Cd(ind0-1)=2/3*Cd(ind0-1);
        Cd(ind0+1)=2/3*Cd(ind0+1);
    end
    clear ind0;
    for j = 1:length(Dd)
        if Dd(j) > 0
            break
        end
    end
    Dd = Dd(j:end);
    tdd = linspace(j,length(Dd)-1,length(Dd));
    if Dd(1) == 0
        %Cd(1) = 2*Cd(2)-Cd(3);
        Dd(1) = (Dd(2))/2;
        Dd(2)= (Dd(2))/2;
        if Dd(2) == 0
            Dd(1) = Dd(3)/3;
            Dd(2) = Dd(3)/3;
            Dd(3) = Dd(3)/3;
        end
    end
    if Dd(end) == 0
        %Cd(1) = 2*Cd(2)-Cd(3);
        temp = 2;
        while Dd(end) <= 0
            for j = length(Dd):-1:length(Dd)-temp+1
                Dd(j) = Dd(end-temp+1)/temp;
            end
            temp = temp+1;
        end
    end
    
    if Dd(1) < 0
        Dd(1) = -Dd(1);
    end
    ind0=find(Dd<=0);
    if isempty(ind0)==0
        Dd(ind0)=1/3*Dd(ind0-1)+1/3*Dd(ind0+1);
        Dd(ind0-1)=2/3*Dd(ind0-1);
        Dd(ind0+1)=2/3*Dd(ind0+1);
    end
    for cc = 1:length(NewTests(1:SHE+21))
        if NewTests(cc) == max(NewTests(1:SHE+21))
            break
        end
    end
    clear ind0;
    
    temp2 = NewTests(1:SHE+21);
    temp2(cc) = [];
    for ccc = 1:length(temp2)
        if temp2(ccc) == max(temp2)
            break
        end
    end
    temp2(ccc) = [];
    


    fprintf('ITERATION %i OF %i\n',its,length(Names))
    tD = linspace(0,length(Cases)-1,length(Cases));
    
    tests = interp1(tD0,avs,tD);
    % [beta, psi, alpha, phi, t, tq, xi, k, A, I0, ag]
    lb = [0, 1, .005, .001,(max(tests(1:SHE))/(0.0050*N))*.7-eps,.4*max(tests(1:SHE)),1,.5,.005*N];
    ub = [6/7.5,6E3,1/7,.05,(max(tests(1:SHE))/(0.0050*N))*.7+eps,.6*max(tests(1:SHE)),800,12,.015*N];
    if TF == 1
        ub(2) = 950;
    end
    if its == 34
        ub(6) = .7*max(tests(1:SHE));
        lb(6) = .5*max(tests(1:SHE));
    end
    %if its == 12
    %    lb(6) = .15*max(NewTests(1:SHE+21));
    %    ub(6) = .35*max(NewTests(1:SHE+21));
    %end
    paramguess = (ub+lb)/2;
    %if its == 53
    %    paramguess = .5*(ub-lb);
    %end
    [paramfit,resnorm] = lsqcurvefit(@objFxn,paramguess,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
    [paramfit,resnorm] = lsqcurvefit(@objFxn,paramfit,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
   % paramfit = [paramfit,21/paramfit(end-1)];
%    lb2 = lb;
 %   ub2 = ub;
  %  lb2(5) = .358 - eps;
   % ub2(5) = .358 + eps;
    %ub2(6) = 2.71E5/sum(StatePops) +eps;
    %lb2(6) = 2.71E5/sum(StatePops) - eps;
    
    %if TF == 1
    %    ub2(2) = 950;
    %end
    %if its == 34
    %    ub2(6) = .7*max(tests(1:SHE));
    %    lb2(6) = .5*max(tests(1:SHE));
    %end
    %paramguess2 = (ub2+lb2)./2;
    %[paramfit2,resnorm2] = lsqcurvefit(@objFxn,paramguess2,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb2,ub2);
    %[paramfit2,resnorm2] = lsqcurvefit(@objFxn,paramfit2,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb2,ub2);
    %paramfit2 = [paramfit2,21/paramfit2(end-1)];
    fits = objFxn(paramfit,tD);
    %Residuals(its,:) = [resnorm,resnorm2];
    %resnorm2-resnorm
    %array2table(Residuals(its,:)) 
    %fits = objFxn(paramfit2,tD);
    Cfit = fits(1:SHE);
    %Cfit2 = fits2(1:SHE);
    Dfit = fits(SHE+1:end);
    %Dfit2 = fits2(SHE+1:end);
    %paramfit = paramfit2;
    k = paramfit(5);
    A = paramfit(6);
    rho = k*((tests/N))./(A/N+((tests/N)))+k0;
    


str = Names(its);
%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
plot(tD(1:SHE+21),Deaths(1:SHE+21),'r*',tD(1:SHE+21),Dfit,'b','linewidth',2)
if TF ~= 1
h=fill([days(StayHomeStart-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeStart-FirstCase)-1],[0,0,max(Deaths(SHE+21),Dfit(SHE+21)),max(Deaths(SHE+21),Dfit(SHE+21))],'k','LineStyle','none');
    h.FaceAlpha=0.1;
end
ylim([0, max(Deaths(SHE+21),Dfit(SHE+21))])
xlim([0,SHE+21])
xticks(round(linspace(0,SHE+21,10)))

xlabel('Days Since First Reported Case up to June 21')
ylabel('Count')
title(str)
legend('Cumulative Death Data','Model Fit','Stay at Home Order in Effect','FontSize',10,'Location','NorthWest')
hold off
set(gca,'FontSize',16)
baseFileName = sprintf('CumDeaths%d',its);
fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure
%hold on
%plot(tD(1:SHE+21),gampdf(tD(1:SHE+21),paramfit(11),paramfit(12)),'b','linewidth',2)
%xlabel('Time Until Death')
%ylabel('density')
%title(str)
%legend('Fit Gamma Distribution PDF')
%xlim([tD(1),tD(SHE+21)])
%hold off
%baseFileName = sprintf('GammaPDF%d',its);
%fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
%saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure
%hold on
%plot(tD(1:SHE+21),gamcdf(tD(1:SHE+21),paramfit(11),paramfit(12)),'b','linewidth',2)
%xlabel('Time Until Death')
%ylabel('cumulative density')
%title(str)
%legend('Fit Gamma Distribution CDF','Location','southeast')
%xlim([tD(1),tD(SHE+21)])
%hold off
%baseFileName = sprintf('GammaCDF%d',its);
%fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
%saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tests=movmean(NewTests,3);
A = paramfit(6);
k = paramfit(5);
tests = interp1(tD0,avs,tD)';
rho = k*((tests/N))./(A/N+((tests/N)))+k0;
% figure
% plot(tD(1:SHE),rho(1:SHE))
%    baseFileName = sprintf('rho%d',its);
%     fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
%     saveas(gca, fullfile(fname, baseFileName), 'png');
% 
 pos = (movmean(Cd,3)./tests(1:end-1))*100;
posFit =(diff(Cfit/lsqwt)'./tests(1:SHE-1))*100;
% for w = 1:length(posFit)
%     if posFit(w) > 100
%         posFit(w) = 100;
%     end
%     if pos(w) > 100
%         pos(w) = 100;
%     end
% end

tP = tD(1:SHE);
%figure
%plot(tP(2:end),pos(1:SHE-1),tP,tests(tP+2))
fig = figure;
   red = [1 0 0];
   blue = [0 0 1];
   left_color = blue;
   right_color = red;
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
   yyaxis left
   hold on
   plot(tP(2:end),pos(1:SHE-1),'b*','linewidth',1)
   plot(tP(2:end),posFit(1:end),'b-','linewidth',2)
   plot(tP,ones(1,length(tP))*5,'b-.')
   if TF ~= 1
   h=fill([days(StayHomeStart-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeStart-FirstCase)-1],[0,0,min(max(pos(1:SHE))+5,100),min(max(pos(1:SHE))+5,100)],'k','LineStyle','none');
    h.FaceAlpha=0.1;
   end
   ylim([0,min(max(pos(1:SHE))+5,100)])
   yticks(round(linspace(0,min(max(pos(1:SHE))+5,100),11),1))
   hold off
   ylabel('Percent Positive Tests')
   yyaxis right
   hold on
   plot(tP,(tests(tP+2)./N)*100,'r--','linewidth',1)   
   plot(tP,(NewTests(tP+2)/N)*100,'r*')
   hold off
   legend({'Positives Inferred From Data','Model Predicted Positives','CDC Positivity Threshold','Stay at Home Order in Effect','New Tests Interpolaton','New Daily Tests Raw Data'},'FontSize',10,'Location','NorthWest')
   ylabel('Percent Population')
   xlabel('Days Since First Reported Case up to May 31')
   xlim([tP(1) tP(end)])
   xticks(round(linspace(tP(1),tP(end),9)))
   yticks(round(linspace(0,1.05*max(((tests(tP+2)./N)*100)),6),2))
   ylim([0,1.05*max(((tests(tP+2)./N)*100))])
   title(str)
   set(gca,'FontSize',16)
    baseFileName = sprintf('Positives%d',its);
    fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
    saveas(gca, fullfile(fname, baseFileName), 'png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure
%hold on
%plot(tD(2:SHE+30),diff(Dfit),'b',tdd,Dd,'r*','linewidth',2)
%if StayHomeEnd ~= NaT
%  h=fill([days(StayHomeStart-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeStart-FirstCase)-1],[0,0,max(Dd(1:SHE+30)),max(Dd(1:SHE+30))],'k','LineStyle','none');
%    h.FaceAlpha=0.1;
%end
%hold off
%xlim([0, SHE+29])
%ylim([0, max(Dd(1:SHE+30))])
%title('Louisiana')
%xlabel('Days Since First Reported Case')
%legend('Model Fit','Inferred Daily Death Totals','Stay at Home Order in Effect')
%ylabel('Count')
%baseFileName = sprintf('DailyDeaths%d',its);
%fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
%saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u0 = [N-paramfit(9)-paramfit(7),paramfit(9),paramfit(7)];
[t,Sol] = ode45(@(t,u)theModel(t,u,paramfit),tD,u0);
infs = Sol(:,3);
R0 = (Sol(:,1)/N)*7.5*paramfit(1);
infsq =  zeros(length(infs),1);
%T = paramfit(9);
z3(1) = paramfit(7)-Cases(1);
for j = 2:length(tD)
    mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1);
    %mid(j-1)=infs(j-1)+infsq(j-1);
    z3(j) = z3(j-1)+(1-rho(j-1))*mid(j-1);
end
figure

hold on
plot(tD(1:SHE),(Cases(1:SHE)/N)*100,'r*',tD(1:SHE),((Cfit/N)*100)/lsqwt,'b',tD(1:SHE),(z3(1:SHE)/N)*100,'b--',tD(1:SHE),(z3(1:SHE)/N)*100+((Cfit/N)*100)/lsqwt,'b-.','linewidth',2)
if TF ~= 1
h=fill([days(StayHomeStart-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeStart-FirstCase)-1],[0,0,(z3(SHE)/N)*100+((Cfit(SHE)/N)*100)/lsqwt,(z3(SHE)/N)*100+((Cfit(SHE)/N)*100)/lsqwt],'k','LineStyle','none');
    h.FaceAlpha=0.1;
end
ylim([0,(z3(SHE)/N)*100+((Cfit(SHE)/N)*100)/lsqwt])
%yticks(round(linspace(0,((z3(SHE)/N)*100+((Cfit(SHE)/N)*100)/lsqwt),11),2))
xlim([tD(1),tD(SHE)])
ylim([0,((z3(SHE)/N)*100+((Cfit(SHE)/N)*100)/lsqwt)])
xticks(round(linspace(tD(1),tD(SHE),9)))
xlabel('Days Since First Reported Case up to May 31')
ylabel('% of population')
title(str)
legend('Reported Cumulative Case Data','Fit','Unreported Case Estimate','True Cumulative Case Estimate','Stay at Home Order in Effect','FontSize',10,'Location','NorthWest')
hold off
set(gca,'FontSize',16)
baseFileName = sprintf('CumCases%d',its);
fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TCT = (z3(1:SHE))+((Cfit))/lsqwt;
%TCTd = diff(TCT);

%figureg
%hold on
%plot(tD(2:SHE),diff(Cfit)/lsqwt,'b',tD(2:SHE),TCTd,'b-.',tdc,Cd,'r*','linewidth',2)
%if TF ~= 1
%  h=fill([days(StayHomeStart-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeStart-FirstCase)-1],[0,0,max(TCTd),max(TCTd)],'k','LineStyle','none');
%    h.FaceAlpha=0.1;
%end
%hold off
%xlim([0, SHE-1])
%ylim([0, max(TCTd)])
%title(str)
%xlabel('Days Since First Reported Case')
%legend('Fit','True Daily Cases Estimate','Inferred Reported Daily Case Totals','Stay at Home Order in Effect')
%ylabel('Count')
%baseFileName = sprintf('DailyCases%d',its);
%fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';a5saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%
paramfit = [paramfit,(TCT(end)/N)*100,(Dfit(end)/N)*100,(Cfit(end)/lsqwt)/TCT(end)];
if its == 7
array2table(paramfit)
end
paramfitsFinal(its,:) = paramfit;
rhos(length(rhos(:,1))-length(rho(1:SHE))+1:end,its) = rho(1:SHE);
R0s(length(rhos(:,1))-length(rho(1:SHE))+1:end,its) = R0(1:SHE);
Cfits(length(Cfits(:,1))-length(Cfit)+1:end,its) = Cfit/lsqwt;
Dfits(length(Dfits(:,1))-length(Dfit)+1:end,its) = Dfit;
URCfits(length(Cfits(:,1))-length(Cfit)+1:end,its) = z3(1:SHE);
Positives(length(Positives(:,1))-length(posFit)+1:end,its) = posFit;
Tests(length(Tests(:,1))-length(tests(1:SHE))+1:end,its) = tests(1:SHE);
%psis = linspace(.01*paramfit(2),paramfit(2),25);
psis = linspace(.1*paramfit(2),paramfit(2),25);
param = paramfit(1:9);
param2 = [paramfit(1:2),0,paramfit(4:9)];
for jj = 1:25
    param(2) = psis(jj);
    param2(2) = psis(jj);
    ff = objFxn(param,tD);
    ff2 = objFxn(param2,tD);
    cf = ff(1:SHE);
    cf2 = ff2(1:SHE);
    u0 = [N-paramfit(9)-paramfit(7),paramfit(9),paramfit(7)];
    [t,Sol] = ode45(@(t,u)theModel(t,u,param),tD,u0);
    [t2,Sol2] = ode45(@(t,u)theModel(t,u,param2),tD,u0);
    infs = Sol(:,3);
    infs2 = Sol2(:,3);
    infsq =  zeros(length(infs),1);
    z4(1) = paramfit(7)-Cases(1);
    z5(1) = z4(1);
for j = 2:length(tD)
    mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1);
    mid2(j-1) =(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs2(j-1:-1:1);
    %mid(j-1)=infs(j-1)+infsq(j-1);
    z4(j) = z4(j-1)+(1-rho(j-1))*mid(j-1);
    z5(j) = z5(j-1)+(1-rho(j-1))*mid2(j-1);
end
    CvaryPsi(jj) = min(100,(100/N)*(((z4(SHE))+((cf(SHE)))/lsqwt)));
    peakC(jj) = (100/N)*max(diff((((z4(1:SHE))+((cf(1:SHE)))/lsqwt))));
    CvaryPsiAlpha(jj) = (100/N)*min(100,100*(((z5(SHE))+((cf2(SHE)))/lsqwt)));
    peakCAlpha(jj) = (100/N)*max(diff(100*(((z5(1:SHE))+((cf2(1:SHE)))/lsqwt))));
    clear infs
end
fprintf('------------------\n')
fprintf('fit alpha:\n')
CvaryPsi(1:5)
fprintf('alpha = 0:\n')
CvaryPsiAlpha(1:5)
fprintf('------------------\n')
PsisWhatIf(:,its) = psis;
CasesWhatif(:,its) = CvaryPsi;
peakCases(:,its) = peakC;
peakCasesAlpha(:,its) = peakCAlpha;
CasesWhatifAlpha(:,its) = CvaryPsiAlpha;
  clear Cases;
    clear Deaths;
    clear NewTests;
    clear SHS;
    clear SHE;
    clear FirstCase;
    clear tW;
    clear tD0;
    clear tD;
    clear tdc;
    clear Cd;
    clear Dd;
    clear tests;
    clear str;
    clear infs;
    clear infsq;
    clear mid;
    clear z3;
    clear Cfit;
    clear Dfit;
    clear temp2;
    clear psis;
    clear cf;
    clear CvaryPsi;
end

save('paramfitsSimplest.mat','paramfitsFinal')
save('rhos.mat','rhos')
save('R0s.mat','R0s')
save('Cfits.mat','Cfits')
save('Dfits.mat','Dfits')
save('URCfits.mat','URCfits')
save('Positives.mat','Positives')
save('Tests.mat','Tests')
save('PsisWhatIf.mat','PsisWhatIf')
save('CasesWhaatIf.mat','CasesWhatif')
save('CasesWhatIfAlpha.mat','CasesWhatifAlpha')
save('peakCases.mat','peakCases')
save('peakCasesAlpha.mat','peakCasesAlpha')
save('Residuals.mat','Residuals')
end