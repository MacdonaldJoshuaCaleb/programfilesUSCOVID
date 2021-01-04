function [] = USStateFits
warning off;
set(0,'DefaultFigureVisible','off')
% import data
paramfits = readtable('paramfits.csv');
paramfits = table2array(paramfits);
CasesByState = readtable('CasesByState.csv');
Dates = table2array(CasesByState(:,1));
CasesByState = table2array(CasesByState(:,2:end));
DeathsByState = readtable('DeathsByState.csv');
DeathsByState = table2array(DeathsByState(:,2:end));
TestingByState = readtable('TestingByState.csv');
TestingByState = table2array(TestingByState(:,2:end));
DatesByState = readtable('ImportantDates.csv');
DatesByState = table2array(DatesByState(:,2:end));
StateNames=["Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","District of Columbia","Florida","Georgia","Guam","Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Northern Mariana Islands","Ohio","Oklahoma","Oregon","Pennsylvania","Puerto Rico","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virgin Islands","Virginia","Washington","West Virginia","Wisconsin","Wyoming"];
Names = StateNames;
StatePops = [4903185,731545,7278717,3017804, 39512223, 5758736, 3565287, 973764, 705749,...
    21477737, 10617423,165768, 1415872, 1787065, 12671821, 6732219, 3155070, 2913314, 4467673,...
    4648794, 1344212, 6045680, 6892503, 9986857, 5639632, 2976149, 6137428, 1068778, 1934408,...
    3080156, 1359711, 8882190, 2096829, 19453561, 10488084, 762062,56882, 11689100, 3956971,...
    4217737,12801989,3193694,1059361, 5148714, 884659, 6829174, 28995881, 3205958, 623989,106977,...
    8535519, 7614893,1792147,5822434, 578759];
 
% remove Louisiana as already handled 
CasesByState(:,20) = [];
DeathsByState(:,20) = [];
TestingByState(:,20) = [];
DatesByState(:,20) = [];
Names(20) = [];
StatePops(20) = [];
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



   function dudt=theModel(t,u,x)
    % u is compartments
    % x is parameter vector, 
    %x = [betag,psig,alphag,phig,Tg,Tqg,xig,kg,Ag]
    beta = x(1);  
    psi = x(2);    
    alpha = x(3); 
    phi = x(4);   
    T = x(5);     
    Tq = x(6); 
    betaq = 0;
    nu = 0;
    phiq = 0;
   
%     tests=interp1(tD0,avs,t);
%     rho = k*tests./(A+tests)+k0;
    dudt = [
        -(beta/N*u(3) + betaq/N*u(4))*u(1)*(1+psi) + alpha*u(2);
        psi*(beta/N*u(3) + betaq/N*u(4))*u(1) - u(2)*(alpha + nu*(beta/N*u(3) + betaq/N*u(4)));
        (beta/N*u(3) + betaq/N*u(4))*((1-phi)*u(1) + nu*(1-phiq)*u(2)) - (1/T)*u(3);
        (beta/N*u(3) + betaq/N*u(4))*(phi*u(1) + nu*phiq*u(2)) - (1/Tq)*u(4)
        ];
end

function z = objFxn(x,tD)
        s0 = N;
        sq0 = 0;
        k = x(8);
        A = x(9);
        z1 = zeros(length(tD),1);
        z2 = zeros(length(tD),1);
        %tests=movmean(NewTests,3);
        tests = interp1(tD0,avs,tD);
        rho = k*(tests./N)./(A/N+tests/N)+k0;
        mdd= zeros(length(tD),1);
        mid= zeros(length(tD),1);
        z1(1)=Cases(1);
        z2(1)=Deaths(1);
%         I0 = Cases(1)/rho(1);
        I0 = x(10);
        ag = x(11);
        bg = 21/ag;
        Iq0 = 0;
     
        T = x(5);     
        Tq = x(6);
        u0 = [s0,sq0,I0,Iq0];
        xi = x(7);
        [t,Sol] = ode45(@(t,u)theModel(t,u,x),tD,u0);
        infs = Sol(:,3);
        infsq =  Sol(:,4);
       % defs = Sol(:,5) + Sol(:,6);
       % z1= Sol(:,5); 
       % mass(1) = (gamcdf(1.5,5.1,.86)-gamcdf(0,5.1,.86))+(gamcdf(1.5,17.8,.45)-gamcdf(0,17.8,.45));
        for j = 2:length(tD)
%             mass(j) =(gamcdf(j+.5,5.1,.86)-gamcdf(j-.5,5.1,.86))+(gamcdf(j+.5,17.8,.45)-gamcdf(j-.5,17.8,.45));
%             defs(j) = mass(j)*defs(j);
            mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1)+(infdist((2:j),1/Tq)-infdist((1:j-1),1/Tq))*infsq(j-1:-1:1);
            mdd(j-1)=(deathdist((2:j),ag,bg)-deathdist((1:j-1),ag,bg))*(infs(j-1:-1:1)+infsq(j-1:-1:1));
            z2(j) = z2(j-1)+mdd(j-1);
            z1(j) = z1(j-1)+rho(j-1)*mid(j-1);
            
        end
%         for j = 1:length(tD)
%             z1(j) = sum(infs(1:j));    
%         end
        z2 = z2*xi;
        z = [lsqwt*z1(1:SHE)',z2(1:SHE+21)'];
end



% Fitting Loop
for its = 1:length(Names)
    Cases = CasesByState(:,its);
    Deaths = DeathsByState(:,its);
    NewTests = TestingByState(:,its);
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
    ind0=find(NewTests<=0);
    if isempty(ind0)==0
        NewTests(ind0)=1/3*NewTests(ind0-1)+1/3*NewTests(ind0+1);
        NewTests(ind0-1)=2/3*NewTests(ind0-1);
        NewTests(ind0+1)=2/3*NewTests(ind0+1);
    end
    
    N = StatePops(its);
    Deaths = Deaths(length(Deaths)-length(Cases)+1:end);
    % get weekly testing averages 
    intv=7;
    Weeks = floor(length(NewTests)/intv);
    avs = zeros(length(Weeks)+1,1);
    for w = 1:Weeks
        avs(w) = mean(NewTests(1+intv*(w-1):intv+intv*(w-1)));
    end
    avs(Weeks+1) = mean(NewTests(Weeks*intv+1:end));
    tW = linspace(0,Weeks,Weeks+1);
    tD0=tW*intv;
    k0=.05;

    % parameter bounds 
    Tg=6;
    betag=3/Tg;
    lsqwt=.05;
    if its == 13
        lsqwt = .1;
    end
    if its == 14
        lsqwt = .1;
    end
    if its == 20
        lsqwt = .01;
    end
    if its == 48
        lsqwt = .01;
    end
    lb = [3/10, 150, .001, .01, 3, 3, .001,.12,min(NewTests(10:SHE))/N,1,.1];
    ub = [6/10,inf,.0995,1.00,10,10,.05,.6,4*max(NewTests)/N,1000,600];

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
    ind0=find(Cd<=0);
    if isempty(ind0)==0
        Cd(ind0)=1/3*Cd(ind0-1)+1/3*Cd(ind0+1);
        Cd(ind0-1)=2/3*Cd(ind0-1);
        Cd(ind0+1)=2/3*Cd(ind0+1);
    end

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

    % [beta, psi, alpha, phi, t, tq, xi, k, A, I0, ag]
    
    paramguess = paramfits(its,1:11);
    %paramguess(9) = paramguess(9)/N;
    
    fprintf('ITERATION %i OF %i\n',its,length(Names))
    tD = linspace(0,length(Cases)-1,length(Cases));
    [paramfit,resnorm] = lsqcurvefit(@objFxn,paramguess,tD(1:SHE+30),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
    %[paramfit,resnorm] = lsqcurvefit(@objFxn,paramfit,tD(1:SHE+30),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+30)'],lb,ub);
    paramfit = [paramfit,21/paramfit(end)];


    fits = objFxn(paramfit,tD);
    Cfit = fits(1:SHE);
    Dfit = fits(SHE+1:end);
    save('fits2.mat','fits')
    k = paramfit(8);
    A = paramfit(9);
    rho = k*(tests./N)./(A/N+tests./N)+k0;



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
xlabel('Days Since First Reported Case')
ylabel('Count')
title(str)
legend('Cumulative Death Data','Model Fit','Stay at Home Order in Effect','Location','southeast')
hold off
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
tests = interp1(tD0,avs,tD)';

pos = (movmean(Cd,3)./tests(2:end))*100;
A = paramfit(9);
k = paramfit(8);
posFit =(diff(Cfit/lsqwt)'./tests(3:SHE+1))*100;
for w = 1:length(posFit)
    if posFit(w) > 100
        posFit(w) = 100;
    end
    if pos(w) > 100
        pos(w) = 100;
    end
end
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
   h=fill([days(StayHomeStart-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeStart-FirstCase)-1],[0,0,100,100],'k','LineStyle','none');
    h.FaceAlpha=0.1;
   end
   hold off
   ylabel('Percent Positive Tests')
   yyaxis right
   hold on
   plot(tP,(tests(tP+2)./N)*100,'r--','linewidth',1)   
   plot(tP,(NewTests(tP+2)/N)*100,'r*')
   hold off
   legend({'Positives Inferred From Data','Model Predicted Positives','CDC Positivity Threshold','Stay at Home Order in Effect','New Tests Interpolaton','New Daily Tests Raw Data'},'Location','northwest')
   ylabel('Percent Population')
   xlabel('Days Since First Reported Case')
   xlim([tP(1) tP(end)])
   title(str)
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

u0 = [N,0,paramfit(10),0];
[t,Sol] = ode45(@(t,u)theModel(t,u,paramfit),tD,u0);
infs = Sol(:,3);
infsq =  Sol(:,4);
T = paramfit(5);
Tq = paramfit(6);
z3(1) = paramfit(10)-Cases(1);
for j = 2:length(tD)
    %mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1)+(infdist((2:j),1/Tq)-infdist((1:j-1),1/Tq))*infsq(j-1:-1:1);
    mid(j-1)=infs(j-1)+infsq(j-1);
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
xlim([tD(1),tD(SHE)])
xlabel('Days Since First Reported Case')
ylabel('% of population')
title(str)
legend('Reported Cumulative Case Data','Fit','Unreported Case Estimate','True Cumulative Case Estimate','Stay at Home Order in Effect','Location','northwest')
hold off
baseFileName = sprintf('CumCases%d',its);
fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TCT = (z3(1:SHE))+((Cfit))/lsqwt;
TCTd = diff(TCT);

figure
hold on
plot(tD(2:SHE),diff(Cfit)/lsqwt,'b',tD(2:SHE),TCTd,'b-.',tdc,Cd,'r*','linewidth',2)
if TF ~= 1
  h=fill([days(StayHomeStart-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeEnd-FirstCase)-1,days(StayHomeStart-FirstCase)-1],[0,0,max(TCTd),max(TCTd)],'k','LineStyle','none');
    h.FaceAlpha=0.1;
end
hold off
xlim([0, SHE-1])
ylim([0, max(TCTd)])
title(str)
xlabel('Days Since First Reported Case')
legend('Fit','True Daily Cases Estimate','Inferred Reported Daily Case Totals','Stay at Home Order in Effect')
ylabel('Count')
baseFileName = sprintf('DailyCases%d',its);
fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%
paramfit = [paramfit,Dfit(end)/TCT(end),resnorm];
paramfitsFinal(its,:) = paramfit;

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
end
save('paramfitsFinal.mat','paramfitsFinal')

  
    
  
end