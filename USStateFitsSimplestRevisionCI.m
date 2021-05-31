function [r] = USStateFitsSimplestRevisionCI(index,sets,remove,noise,option) 
warning off;
if index == 1
close all
weight = 5.91;
end
set(0,'DefaultFigureVisible','on')
% import data
params = readtable('DeathsAsChinaParams.csv');
params = table2array(params(:,2:9));
params(:,1) = params(:,1)./7.5;
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
paramfits = zeros(55,16);
cumcases = zeros(131,55);
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
%FirstCases(52)

function [d1]=getData(CFit,nLevel,sets)
    d1 = zeros(sets,length(CFit));
    for k = 1:sets
        for j = 1:length(CFit)
            d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
            while d1(k,j) < 0
                d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
            end
            %d1(k,j) = ceil(d1(k,j));
        end
    end
end

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
        (beta/N*u(3))*u(1) - (1/T)*u(3); (beta/N)*u(3)*u(1) 
        ];
   end
function z = objFxn(x,tD)
        s0 = N;
        %sq0 = x(9);
        k = x(5);
        A = x(6);
        beta = x(1);
        psi = x(2);
        z1 = zeros(length(tD),1);
        z2 = zeros(length(tD),1);
        %tests=movmean(NewTests,3);
        tests = interp1(tD0,avs,tD);
        %tests = ceil(tests);
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
        %sq0 = x(9);
        %T = x(9);     
        
        
         u0 = [s0-I0,0,I0,(beta/N)*I0*(N-I0)];
        xi = x(4);
        [t,Sol] = ode45(@(t,u)theModel(t,u,x),tD,u0);
        infs = Sol(:,3);
        infsq =  zeros(length(infs),1);
         tmp = zeros(length(tD),1);
                  tmp(1) = Sol(1,end);
         tmp(2:end) = Sol(2:end,end)-Sol(1:end-1,end);

       % defs = Sol(:,5) + Sol(:,6);
       % z1= Sol(:,5); 
       % mass(1) = (gamcdf(1.5,5.1,.86)-gamcdf(0,5.1,.86))+(gamcdf(1.5,17.8,.45)-gamcdf(0,17.8,.45));
        for j = 2:length(tD)
%             mass(j) =(gamcdf(j+.5,5.1,.86)-gamcdf(j-.5,5.1,.86))+(gamcdf(j+.5,17.8,.45)-gamcdf(j-.5,17.8,.45));
%             defs(j) = mass(j)*defs(j);
            mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*tmp(j-1:-1:1);
            mdd(j-1)=(deathdist((2:j),ag,bg)-deathdist((1:j-1),ag,bg))*(tmp(j-1:-1:1));
            z2(j) = z2(j-1)+mdd(j-1);
            z1(j) = z1(j-1)+rho(j-1)*mid(j-1);
            % get rid z2(j-1), z1(j-1)
        end
%         for j = 1:length(tD)
%             z1(j) = sum(infs(1:j));    
%         end
        z2 = z2*xi;
        cs = z1;
        ds = z2;
        cf = diff(cs(1:SHE));
        %cf(1) = Cd(1);
        df = diff(ds(1:SHE+21));
        %df(1) = Dd(1);
        %(diff(cs(1:SHE))'./tests(1:SHE-1))*lsqwt3
        z = [lsqwt1*cs(1:SHE)',lsqwt2*ds(1:SHE+21)'];
end

% Fitting Loop
its = index;
 %   set(0,'DefaultFigureVisible','off')
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
     Cd = diff(Cases);
     WeeksC = floor(length(Cd)/intv);
    Dd = diff(Deaths);
    avs = zeros(length(Weeks)+1,1);
    avsC = zeros(length(Weeks)+1,1);
    for w = 1:Weeks
        avs(w) = mean(NewTests(1+intv*(w-1):intv+intv*(w-1)));
    end
    for w = 1:WeeksC
         avsC(w) = mean(Cd(1+intv*(w-1):intv+intv*(w-1)));
    end
    avs(Weeks+1) = mean(NewTests(Weeks*intv+1:end));
    avsC(WeeksC+1) = mean(Cd(WeeksC*intv+1:end));
    tW = linspace(0,Weeks,Weeks+1);
    tWC = linspace(0,WeeksC,WeeksC+1);
    tD0=tW*intv;
    tD0C = tWC*intv;
    k0=.03;
    if index == 2
        k0 = .07;
    end
    if index == 36
        k0 = .05;
    end
    

   
  

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
    
   % Dd = Dd(1:SHE+20);
   % Cd = Cd(1:SHE-1);
    lsqwt1 = 1/max(Cases(1:SHE));
    lsqwt2 = 2.5/max(Deaths(1:SHE+21));
    fprintf('ITERATION %i OF %i\n',its,length(Names))
    tD = linspace(0,length(Cases)-1,length(Cases));
    
    tests = interp1(tD0,avs,tD);
    %tests = ceil(tests);
    I0g = (Cases(1)/.005)*T;
    if index == 6
        I0g = (Cases(1)/.03)*T;
    end
    if index == 14
        I0g = (Cases(1)./.03)*T;
    end
    if index == 36
        I0g = 2*(Cases(1)./.03)*T;
    end
    if index == 17
        I0g = (Cases(1)./.03)*T;
    end
    if index == 26
        I0g = (Cases(1)./.03)*T;
    end
    % [beta, psi, alpha, phi, t, tq, xi, k, A, I0, ag]
    lb = [0,     1,  1/200,  .001,       min(1-eps,(max(tests(1:SHE))/(0.0050*N))*1.1-eps),.4*max(tests(1:SHE)),1,.5];
    ub = [6/7.5,1000, 1/7,.05, min(1,(max(tests(1:SHE))/(0.0050*N))*1.1+eps),.6*max(tests(1:SHE)),I0g,25];
    if its == 17
        array2table([lb,ub])
    end
    if its == 37
        %lb(5) = lb(5);
        %ub(5) = ub(5);
        lb(6) = .8*max(tests(1:SHE));
        ub(6) = max(tests(1:SHE));
        ub(7) = 5;
        lb(5) = .9;
        ub(5) = 1;
    end
    if its == 12
        ub(7) = 5;
    end
%     if its == 34
%         lb(6) = .4*max(tests(1:SHE));
%         ub(6) = .6*max(tests(1:SHE));
%     end
%     if its == 37
%         ub(7) = 200;
%         ub(2) = 4000;
%         lb(5) = (max(tests(1:SHE))/(0.0050*N))*.9-eps;
%         ub(5) = (max(tests(1:SHE))/(0.0050*N))*.9-eps;
%     end
%    
    %
    %if its == 12
    %    lb(6) = .15*max(NewTests(1:SHE+21));
    %    ub(6) = .35*max(NewTests(1:SHE+21));
    %end
    paramguess = params(its,:);
    if its == 37
        paramguess = (lb+ub)./2;
    end
    if its == 6
        paramguess(7) = 100;
        array2table(paramguess)
    end
    if its == 14
        paramguess(7) = 200;
        array2table(paramguess)
    end
    if its == 17
        paramguess(7) = 25;
        paramguess(7) = 1/8;
        array2table(paramguess)
        array2table([lb;ub])
    end
    if its == 26
        paramguess(7) = 200;
        paramguess(3) = 1/10;
    end
     if its == 36
%         paramguess(7) = 500;
      paramguess(3) = 1/8;
%         array2table(paramguess)
    end
    %end
   p = movmean(diff(Cases),3)./tests(1:end-1)';

    pos = p(1:SHE-1);
    for jj = 1:length(pos)
        if pos(jj) > 1
            pos(jj) = 1;
        end
    end
        if its == 34
            lsqwt2 = 1.25./Deaths(SHE+21);
        end
         if its == 17
             lsqwt2 = .75./Deaths(SHE+21);
         end
         if its == 26
              lsqwt2 = 1./Deaths(SHE+21);
         end
         if its == 36
              lsqwt2 = 5.91/Deaths(SHE+21);
         end
%         if its == 36
%             paramguess = [0.1957       58.842        1/10      0.012036      0.70941      710.17       121.68       5.2671 ];
%         end
    [paramfitguess,resnorm] = lsqcurvefit(@objFxn,paramguess,tD(1:SHE+21),[lsqwt1*Cases(1:SHE)',lsqwt2*Deaths(1:SHE+21)'],lb,ub);
    lb(6) = .85*paramfitguess(6);
    ub(6) = 1.15*paramfitguess(6);
    if its ~= 17
    if its ~= 14
    if its ~= 6
    lb(7) = .5*paramfitguess(7);
    ub(7) = 2*paramfitguess(7);
    end
    end
    end
    if its == 34
        lb(6) = .85*lb(6);
        ub(6) = 1.15*ub(6);
    end
    if its == 32
        lb(6) = .85*lb(6);
        ub(6) = 1.15*ub(6);
    end

%r=paramfit;
if its == 45
    paramfitguess(3) = .1181;
end
    [paramfit,resnorm] = lsqcurvefit(@objFxn,paramfitguess,tD(1:SHE+21),[lsqwt1*Cases(1:SHE)',lsqwt2*Deaths(1:SHE+21)'],lb,ub);
     ub(6) = 1.15*paramfit(6);
   lb(6) = .85*paramfit(6);
    if its == 34
        lb(6) = .85*lb(6);
        ub(6) = 1.15*ub(6);
    end
      if its == 32
        lb(6) = .85*lb(6);
        ub(6) = 1.15*ub(6);
    end
   lb(7) = .5*paramfit(7);
        ub(7) = 2*paramfit(7);
        if lb(7) < 1
            lb(7) = 1;
            ub(7) = (Cases(1)/.03)*T;
            array2table(ub)
        end
check = [6,18,20,24];
for c = 1:length(check)
    if its == check(c)
        lb(7) = 1;
        ub(7) = 2000;
        paramfit(7) = 500;
        break
    end
end
 check2 = [6,12,17,37,40,42];
 for c = length(check2)
 if its == check2(c)
     lb(7) = 1;
 end
 end
% 
     if its == 38
         lb(7) = 1;
         ub(7) = 3400;
         paramfit = (lb+ub)./2;
         paramfit(7) = 1500;
         
     end
         if its == 41
         ub(7) = 2500;
         lb(7) = 1;
         paramfit = (lb+ub)./2;
         paramfit(7) = 1500;
         
         end
     if ub(7) < lb(7)
         lb(7) = max(.25*ub(7),Cases(1));
     end
%array2table([lb;ub])
array2table([lb;ub])
   [paramfit,resnorm] = lsqcurvefit(@objFxn,paramfit,tD(1:SHE+21),[lsqwt1*Cases(1:SHE)',lsqwt2*Deaths(1:SHE+21)'],lb,ub);
        
    fits = objFxn(paramfit,tD(1:SHE+21));
    cfit = fits(1:SHE)./lsqwt1;
    dfit = fits(SHE+1:end)./lsqwt2;
    %pfit = fits(2*SHE+22:end)./lsqwt3;
    
    
    figure
    hold on
    plot(0:length(Deaths(1:SHE+21))-1,(Deaths(1:SHE+21)./N)*100,'b+','linewidth',2)
    plot(0:length(Deaths(1:SHE+21))-1,(dfit./N)*100,'b','linewidth',2)
    %plot(1:length(Deaths(1:SHE+21))-1,(movmean(Dd(1:SHE+20),3)./N)*100,'k')
    hold off
    ylabel('Cum. Deaths (% population)')
    xlabel('Days Since First Reported Case to June 21')
    xlim([0,SHE+20])
    %legend('Reported Deaths', 'Deaths Fit','Fontsize',12,'Location','NorthWest')
     set(gca,'FontSize',16)
%     baseFileName = sprintf('CumDeaths%d',its);
% fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
% saveas(gca, fullfile(fname, baseFileName), 'png');
       I0 = paramfit(7);
       
        %T = x(9);     
        
       
    psi = paramfit(2);
    beta = paramfit(1);
    I0 = paramfit(7);
     u0 = [N-I0,0,I0,(beta/N)*I0*(N-I0)];
    [t,Sol] = ode45(@(t,u)theModel(t,u,paramfit),tD,u0);
%     figure
%    plot(0:1:SHE-1,cfit./Sol(1:SHE,4)','b','linewidth',2)
   
    txt = 'Cum. Case Est.';
    txt2 = 'Daily Case Est.';
   % txt3 = 'Rep. Case Fit';
       [M I] = max(Sol(1:SHE,3));
     
   fig = figure;
   red = [1 0 0];
   blue = [0 0 1];
   left_color = blue;
   right_color = red;
   set(fig,'defaultAxesColorOrder',[left_color; right_color]);
   yyaxis left
   hold on
    plot(0:length(Cases(1:SHE))-1,(Sol(1:length(Cases(1:SHE)),end)./N)*100,'b-.','linewidth',2)
   % text(I,(Sol(I,4)./N)*100,txt,'Fontsize',12)
    plot(0:length(Cases(1:SHE))-1,(Cases(1:SHE)./N)*100,'k+','linewidth',2)
    plot(0:length(Cases(1:SHE))-1,(cfit./N)*100,'b-','linewidth',2)
  %  text(15,1.1*(cfit(60)./N)*100,txt3,'Fontsize',12)
    plot(t(1:SHE),(Sol(1:SHE,3)./N)*100,'b:','linewidth',2)
  %  text(I+5,1.2*(Sol(I+5,3)./N)*100,txt2,'Fontsize',12)
   hold off
   ylabel('Cases (% population)')
      yyaxis right
   hold on
    plot(t(1:SHE),(Sol(1:SHE,2)./N)*100,'r--','linewidth',2)
   hold off
   ylabel('Self-Quar. (% population)')
   xlabel('Days Since First Reported Case to May 31')
   xlim([0,SHE-1])
   set(gca,'FontSize',16)
  % legend({'Cum. Case Estimate','Reported Cases','Reported Cases Fit','Daily Cases'},'Fontsize',16,'Location','EastOutside')
% baseFileName = sprintf('CumCases%d',its);
% fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
% saveas(gca, fullfile(fname, baseFileName), 'png');

posData = movmean(Cd(1:SHE-1),3)'./tests(1:SHE-1);
posFit = (diff(cfit)./tests(1:SHE-1));
for jj = 1:length(posData)
    if posData(jj) > 1
        posData(jj) = 1;
    end
    if posFit(jj) > 1
        posFit(jj) = 1;
    end
end


 
    

    R0 = (Sol(:,1)/N)*7.5*paramfit(1);
     k = paramfit(5);
        A = paramfit(6);

        rho = k*((tests/N))./(A/N+((tests/N)))+k0;
        
        
  fig2 = figure;
   red = [1 0 0];
   blue = [0 0 1];
   left_color = blue;
   right_color = red;
   set(fig2,'defaultAxesColorOrder',[left_color; right_color]);
   yyaxis left
   hold on
  plot(1:SHE-1,posFit*100,'b-','linewidth',2)
  plot(1:SHE-1,posData*100,'b+','linewidth',2)
   hold off
   ylabel('positive tests (%)')
      yyaxis right
   hold on
      plot(1:SHE-1,rho(1:SHE-1),'r-','linewidth',2)
    %plot(1:SHE-1,(NewTests(1:SHE-1)/N)*100,'r*','linewidth',2)
   hold off
   ylabel('Daily Tests (% population)')
   xlabel('Days Since First Reported Case to May 31')
   xlim([0,SHE-1])
   set(gca,'FontSize',16)
   title(Names(its))
  set(gca,'FontSize',16)
%     baseFileName = sprintf('Positives%d',its);
%     fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
%     saveas(gca, fullfile(fname, baseFileName), 'png');

PsisWhatIf = zeros(25,55);
CasesWhatif = zeros(25,55);
CasesWhatifAlpha = zeros(25,55);


psis = linspace(.1*paramfit(2),paramfit(2),25);
 param = paramfit(1:8);
 param2 = [paramfit(1:2),0,paramfit(4:8)];
 u0Alt = [N-param(7),0,param(7),(param(1)/N)*param(7)*(N-param(7))];
 for jj = 1:25
    param(2) = psis(jj);
     param2(2) = psis(jj);
%     ff = objFxn(param,tD);
%     ff2 = objFxn(param2,tD);
%     cf = ff(1:SHE);
%     cf2 = ff2(1:SHE);
 %    u0 = [N-paramfit(9)-paramfit(7),paramfit(9),paramfit(7)];
     [t,Sol2] = ode45(@(t,u)theModel(t,u,param),tD,u0Alt);
     [t2,Sol3] = ode45(@(t,u)theModel(t,u,param2),tD,u0Alt);
     CvaryPsi(jj) = (100/N)*(Sol2(SHE,end));
     peakC(jj) = (100/N)*max(Sol2(1:SHE,3));
     CvaryPsiAlpha(jj) = (100/N)*(Sol3(SHE,end));
     peakCAlpha(jj) = (100/N)*max(Sol3(1:SHE,3));
 end


if option == 0
r = [psis',CvaryPsi',peakC',CvaryPsiAlpha',peakCAlpha'];
end





paramfits(its,:) = [paramfit, (Sol(length(Cases(1:SHE)),end)./N)*100,  Sol(length(Cases(1:SHE)),end)./cfit(end), (M/N)*100,I,-days(FirstCases(52)-FirstCase),(Sol(SHE,3)./N)*100,(dfit(end)./N)*100,(cfit(end)./N)*100];
paramfit = [paramfit,(Sol(length(Cases(1:SHE)),end)./N)*100,Sol(length(Cases(1:SHE)),end)./cfit(end),(M/N)*100,I,0,(Sol(SHE,3)./N)*100,(dfit(end)./N)*100,(cfit(end)./N)*100];
 %array2table(paramfit-paramfits(its,:))
% if option == 0
%     r = paramfit;


Sol1 = Sol;
%  for l = 1:200
%      Cstar(1) = poissrnd(cfit(1));
%      Dstar(1) = poissrnd(dfit(1));
%      for jj = 2:length(cfit)
%          Cstar(jj) = Cstar(jj-1)+poissrnd((cfit(jj)-cfit(jj-1)));
%      end
%      for jj = 2:length(dfit)
%          Dstar(jj) = Dstar(jj-1)+poissrnd(dfit(jj)-dfit(jj-1));
%      end
%      CSstar(l,:) = Cstar;
%      DSstar(l,:) = Dstar;
%  end
%  
%  figure 
%  hold on
%  plot(0:1:130,(cfit),'k','linewidth',2)
%  plot(0:1:130,(CSstar'),'b')
%  plot(0:1:130,(Cases(1:131)),'r*')
%  hold off
%  
%  figure 
%  hold on
%  plot(0:1:130+21,(dfit),'k','linewidth',2)
%  plot(0:1:130+21,(DSstar'),'b')
%  plot(0:1:130+21,(Deaths(1:SHE+21)),'r*')
%  hold off
 
 

%sets = 200;
if option == 1
  
SynthCumC = zeros(sets,SHE);
SynthCumD = zeros(sets,SHE+21);
for jj = 1:sets
    if mod(jj,.1*sets) == 0
        fprintf('-----------------------------------------------------------\n')
        fprintf('WORKING ON SET %i OF % i\n',jj,sets)
        fprintf('-----------------------------------------------------------\n')
    end
Cs = getData(diff(cfit),noise,1);
Ds = getData(diff(dfit),noise,1);
    for kk = 2:SHE+21
        if kk == 2
            SynthCumC(jj,kk-1) = Cases(1);
            SynthCumD(jj,kk-1) = Deaths(1);
        end
        if kk <= SHE
        SynthCumC(jj,kk) = sum(Cs(1:kk-1));
        end
        SynthCumD(jj,kk) = sum(Ds(1:kk-1));
    end
    
% SynthCumC(jj,1) = poissrnd(cfit(1));
% SynthCumD(jj,1) = dfit(1);
% for kk=2:SHE+21
%     if kk <= SHE
%     SynthCumC(jj,kk) = SynthCumC(jj,kk-1) + poissrnd(cfit(kk)-cfit(kk-1));
%     end
%     SynthCumD(jj,kk) = SynthCumD(jj,kk-1) + poissrnd(dfit(kk)-dfit(kk-1));
% end
%SynthCumD = SynthCumD*paramfit(4);
%  lsqwt1 = 1./SynthCumC(jj,end);
%  lsqwt2 = 2.5./SynthCumD(jj,end);
%lb(6) = .85*paramfit(6);
%ub(6) = 1.15*paramfit(6);
%lb(5) = paramfit(5)-.02;
%ub(5) = paramfit(5)+.02;
% lb(7) = .25*paramfit(7);
% ub(7) = 4*paramfit(7);

lsqwt1 = 1/SynthCumC(jj,end);
lsqwt2 = 2.5/SynthCumD(jj,end);
if its == 34
lsqwt2 = 1.25/SynthCumD(jj,end);
end
if its == 17
    lsqwt1 = 1/SynthCumC(jj,end);
    lsqwt2 = .75/SynthCumD(jj,end);
end
if its == 26
    lsqwt1 = 1/SynthCumC(jj,end);
    lsqwt2 = 1/SynthCumD(jj,end);
end
if its == 36
    lsqwt1 = 1/SynthCumC(jj,end);
    lsqwt2 = 5.91/SynthCumD(jj,end);
end
lb(7) = .5*paramfit(7);
ub(7) = 1.5*paramfit(7);
% if its == 6
%     ub(7) = 1200;
% end

%if its ~= 36
lb(6) = .85*paramfit(6);
ub(6) = 1.15*paramfit(6);
%end
%if its ~= 17
lb(5) = paramfit(5)-.02;
ub(5) = paramfit(5)+.02;
%end
if jj == 1
array2table([lb;ub])
end
    [paramfitInt,resnorm] = lsqcurvefit(@objFxn,paramfit(1:8),tD(1:SHE+21),[lsqwt1*SynthCumC(jj,:),lsqwt2*SynthCumD(jj,:)],lb,ub);
    fits = objFxn(paramfitInt,tD(1:SHE+21));
    cfits = fits(1:SHE)./lsqwt1;
    dfits = fits(SHE+1:end)./lsqwt2;
    I0 = paramfitInt(7);
    psi = paramfitInt(2);
    beta = paramfitInt(1);
     u0 = [N-I0,0,I0,(beta/N)*I0*(N-I0)];
    [t,Sol] = ode45(@(t,u)theModel(t,u,paramfitInt),tD,u0);
    [M I] = max(Sol(1:SHE,3));
    paramfitInt = [paramfitInt,(Sol(length(Cases(1:SHE)),end)./N)*100,Sol(length(Cases(1:SHE)),end)./cfits(end),(M/N)*100,I,0,(Sol(SHE,3)./N)*100,(dfits(end)./N)*100,(cfits(end)./N)*100];
     R0s(jj,:) = (Sol(1:SHE,1)/N)*7.5*paramfitInt(1);
     k = paramfitInt(5);
     A = paramfitInt(6);
     rhoInt = k*((tests/N))./(A/N+((tests/N)))+k0;
     rhos(jj,:) = rhoInt(1:SHE);
     ccfits(jj,:) = cfits;
     ddfits(jj,:) = dfits;
     CumCfits(jj,:) = Sol(1:SHE,4);
     Quars(jj,:) = Sol(1:SHE,2);
     DCs(jj,:) = Sol(1:SHE,3);
     dcc = diff(cfits');
     dcc = dcc';
     posFit = (dcc./tests(1:SHE-1));
     for tt = 1:length(posFit)
         if posFit(tt) > 1
             posFit(tt) = 1;
         end
     end
     posFits(jj,:) = posFit;
     paramfitInts(jj,:) = paramfitInt; 
     REs(jj,:) = abs(paramfitInt-paramfit)./paramfit;
    %array2table(paramfitInt)
end

% figure 
% hold on
% plot(0:1:SHE-1,100*SynthCumC./N,'b')
% plot(0:1:SHE-1,100*Cases(1:SHE)./N,'r*')
% hold off
% title('Generated Datasets, Cum. Cases')
% xlabel('Days Since First Reported Case to May 31')
% ylabel('% pop.')
% set(gca,'FontSize',16)
% 
% size(diff(SynthCumC'))
% length(1:1:SHE)
% figure 
% hold on
% plot(1:1:SHE-1,100*diff(SynthCumC')./N,'b')
% plot(1:1:SHE-1,100*diff(Cases(1:SHE))./N,'r*')
% hold off
% title('Generated Datasets, Daily Cases')
% xlabel('Days Since First Reported Case to May 31')
% ylabel('% pop.')
% set(gca,'FontSize',16)
% 
% figure 
% hold on
% plot(0:1:SHE-1+21,100*SynthCumD./N,'b')
% plot(0:1:SHE-1+21,100*Deaths(1:SHE+21)./N,'r*')
% hold off
% title('Generated Datasets, Cum. Deaths')
% xlabel('Days Since First Reported Case to May 31')
% ylabel('% pop.')
% set(gca,'FontSize',16)
% 
% figure 
% hold on
% plot(1:1:SHE-1+21,100*diff(SynthCumD')./N,'b')
% plot(1:1:SHE-1+21,100*diff(Deaths(1:SHE+21))./N,'r*')
% hold off
% title('Generated Datasets, Daily Deaths')
% xlabel('Days Since First Reported Case to June 21')
% ylabel('% pop.')
% set(gca,'FontSize',16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1 = [0:1:SHE-1+21,fliplr(0:1:SHE-1+21)];
x2 = [0:1:SHE-1,fliplr(0:1:SHE-1)];
S = 7;
if index == 46
    S = 8;
end
x3 = [S-1:SHE-1,fliplr(S-1:SHE-1)];
x4 = [S:1:SHE-1,fliplr(S:1:SHE-1)];
%remove = 5;
srhos = sort(rhos);
rhoLB = srhos(1+remove,:);
rhoUB = srhos(end-remove,:);

% figure
% hold on
% inBetween = [rhoLB,fliplr(rhoUB)];
% h = fill(x2, inBetween,'b','LineStyle','none');
% set(h,'facealpha',.4)
% %plot(0:1:130,rhoLB,'b')
% %plot(0:1:130,rhoUB,'b')
% plot(0:1:SHE-1,rho(1:SHE),'b','linewidth',2)
% xlim([0 SHE-1])
% hold off
% %title('Generated Datasets, Cum. Cases')
% ylabel('Ascertainment Rate')
% xlabel('Days Since First Reported Case to May 31')
% set(gca,'FontSize',16)

sR0s = sort(R0s);
R0LB = sR0s(1+remove,:);
R0UB = sR0s(end-remove,:);

% figure
% hold on
% inBetween = [R0LB,fliplr(R0UB)];
% h = fill(x2, inBetween,'b','LineStyle','none');
% set(h,'facealpha',.4)
% plot(0:1:SHE-1,R0(1:SHE),'b','linewidth',2)
% xlim([0 SHE-1])
% hold off
% ylabel('Time Variable R_e')
% xlabel('Days Since First Reported Case to May 31')
% set(gca,'FontSize',16)


sCC = sort(CumCfits);
CCLB = 100*sCC(1+remove,:)./N;
CCUB = 100*sCC(end-remove,:)./N;
sRC = sort(ccfits);
RCLB = 100*sRC(1+remove,:)./N;
RCUB = 100*sRC(end-remove,:)./N;
sQ = sort(Quars);
QLB = 100*sQ(1+remove,:)./N;
QUB = 100*sQ(end-remove,:)/N;
sDC = sort(DCs);
DCLB = 100*sDC(1+remove,:)./N;
DCUB = 100*sDC(end-remove,:)./N;

 Casefig = figure;
   red = [1 0 0];
   blue = [0 0 1];
   left_color = blue;
   right_color = red;
 set(Casefig,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
hold on
inBetween = [CCLB,fliplr(CCUB)];
h = fill(x2, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
plot(0:1:SHE-1,100*Sol1(1:SHE,4)./N,'b-.','linewidth',2)
%plot(0:1:130,(CumCfits./N)*100,'-.','linewidth',2)
inBetween = [RCLB,fliplr(RCUB)];
h = fill(x2, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
plot(0:1:SHE-1,100*cfit./N,'b-','linewidth',2)
%plot(0:1:130,(ccfits./N)*100,'-','linewidth',2)
plot(0:1:SHE-1,(Cases(1:SHE)./N)*100,'k+')
%plot(0:1:130,(DCs./N)*100,':','linewidth',2)
inBetween = [DCLB,fliplr(DCUB)];
h = fill(x2, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
plot(0:1:SHE-1,100*Sol1(1:SHE,3)./N,'b:','linewidth',2)
hold off
ylabel('Cases (% pop.)')
yyaxis right
hold on
inBetween = [QLB,fliplr(QUB)];
h = fill(x2, inBetween,'r','LineStyle','none');
set(h,'facealpha',.4)
plot(0:1:SHE-1,100*Sol1(1:SHE,2)./N,'r--','linewidth',2)
%plot(0:1:130,(Quars./N)*100,'-','linewidth',2)
hold off
ylabel('Self-Quar. (% pop.)')
xlabel('Days Since First Reported Case to May 31')
xlim([0, SHE-1])
set(gca,'FontSize',16)
baseFileName = sprintf('CumCases%d',its);
fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
saveas(gca, fullfile(fname, baseFileName), 'png');
sD = sort(ddfits);
DLB = 100*sD(1+remove,:)./N;
DUB = 100*sD(end-remove,:)./N;

figure
hold on
inBetween = [DLB,fliplr(DUB)];
h = fill(x1, inBetween,'k','LineStyle','none');
set(h,'facealpha',.4)
%plot(0:1:130+21,(ddfits./N)*100,'b','linewidth',2)
plot(0:1:SHE-1+21,(Deaths(1:SHE+21)./N)*100,'k+')
plot(0:1:SHE-1+21,100*dfit./N,'k','linewidth',2)
xlim([0 SHE-1+21])
hold off
ylabel('Deaths (% pop.)')
xlabel('Days Since First Reported Case to June 21')
title(Names(its))
set(gca,'FontSize',16)
baseFileName = sprintf('CumDeaths%d',its);
fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
saveas(gca, fullfile(fname, baseFileName), 'png');
%for jj = 1:130
%    if posFit(jj) > 1
%        posFit(jj) = 1;
%    end
%end

sP = sort(posFits);
PLB = 100*sP(1+remove,S:end);
PUB = 100*sP(end-remove,S:end);
srhos = sort(rhos);
rhoLB = srhos(1+remove,S:end);
rhoUB = srhos(end-remove,S:end);
 Posfig = figure;
   red = [1 0 0];
   blue = [0 0 1];
   left_color = blue;
   right_color = red;
 set(Posfig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
hold on
plot(S:SHE-1,posData(S:end)*100,'b+','linewidth',1)
inBetween = [PLB,fliplr(PUB)];
h = fill(x4, inBetween,'b','LineStyle','none');
set(h,'facealpha',.4)
plot(S:SHE-1,posFit(S:SHE-1)*100,'b-','linewidth',2)
ylabel('postive tests (%)')
hold off

yyaxis right
hold on
inBetween = [rhoLB,fliplr(rhoUB)];
h = fill(x3, inBetween,'r','LineStyle','none');
set(h,'facealpha',.4)
%plot(0:1:130,rhoLB,'b')
%plot(0:1:130,rhoUB,'b')
plot(S-1:1:SHE-1,rho(S:SHE),'r-','linewidth',2)
xlim([S-1 SHE-1])
hold off
%title('Generated Datasets, Cum. Cases')
ylabel('Ascertainment Rate')
xlabel('Days Since First Reported Case to May 31')
xlim([S, SHE-1])
set(gca,'FontSize',16)
     baseFileName = sprintf('Positives%d',its);
     fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
     saveas(gca, fullfile(fname, baseFileName), 'png');


DRC = diff(ccfits');
DRC = DRC';
sDRC = sort(DRC);
DRCLB = 100*sDRC(1+remove,:)./N;
DRCUB = 100*sDRC(end-remove,:)./N;

% figure
% hold on
% inBetween = [DRCLB,fliplr(DRCUB)];
% h = fill(x3, inBetween,'k','LineStyle','none');
% set(h,'facealpha',.4)
% plot(1:1:SHE-1,100*diff(cfit)./N,'b','linewidth',2)
% plot(1:1:SHE-1,100*Cd(1:SHE-1)./N,'k+')
% inBetween = [DCLB,fliplr(DCUB)];
% h = fill(x2, inBetween,'b','LineStyle','none');
% set(h,'facealpha',.4)
% plot(0:1:SHE-1,100*Sol1(1:SHE,3)./N,'k:','linewidth',2)
% hold off
% ylabel('Daily Cases (% pop.)')
% xlabel('Days Since First Reported Case to May 31')
% xlim([0 SHE-1])
% set(gca,'FontSize',16)


CA = ccfits./CumCfits;
sCA = sort(CA);
CALB = sCA(1+remove,1:end);
CAUB = sCA(end-remove,1:end);



Ascfig = figure;
%    red = [1 0 0];
%    blue = [0 0 1];
%    left_color = blue;
%    right_color = red;
%  set(Ascfig,'defaultAxesColorOrder',[left_color; right_color]);
hold on
inBetween = [CALB,fliplr(CAUB)];

h = fill(x2, inBetween,'r','LineStyle','none');
set(h,'facealpha',.4)
%plot(0:1:SHE-1,ccfits./CumCfits,'b')
plot(0:1:SHE-1,cfit(1:end)./Sol1(1:SHE,4)','r-','linewidth',2)
hold off
ylabel('Cum. Ascertainment Ratio')
xlabel('Days Since First Reported Case to May 31')
%xlabel('Days Since First Reported Case to May 31')
xlim([0 SHE-1])
set(gca,'FontSize',16)
    baseFileName = sprintf('Ascers%d',its);
     fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
     saveas(gca, fullfile(fname, baseFileName), 'png');




%  Ascfig = figure;
%    red = [1 0 0];
%    blue = [0 0 1];
%    left_color = blue;
%    right_color = red;
%  set(Ascfig,'defaultAxesColorOrder',[left_color; right_color]);
% yyaxis left
% hold on
% inBetween = [rhoLB,fliplr(rhoUB)];
% h = fill(x2, inBetween,'b','LineStyle','none');
% set(h,'facealpha',.4)
% %plot(0:1:130,rhoLB,'b')
% %plot(0:1:130,rhoUB,'b')
% plot(0:1:SHE-1,rho(1:SHE),'b-','linewidth',2)
% xlim([0 SHE-1])
% hold off
% %title('Generated Datasets, Cum. Cases')
% ylabel('Ascertainment Rate')
% xlabel('Days Since First Reported Case to May 31')
% set(gca,'FontSize',16)
% 
% yyaxis right
% hold on
% inBetween = [CALB,fliplr(CAUB)];
% size(inBetween)
% size(x4)
% h = fill(x4, inBetween,'r','LineStyle','none');
% set(h,'facealpha',.4)
% %plot(0:1:SHE-1,ccfits./CumCfits,'b')
% plot(4:1:SHE-1,cfit(5:end)./Sol1(5:SHE,4)','r-','linewidth',2)
% hold off
% ylabel('Cum. Ascertainment Ratio')
% %xlabel('Days Since First Reported Case to May 31')
% xlim([0 SHE-1])
% set(gca,'FontSize',16)
%     baseFileName = sprintf('Ascers%d',its);
%      fname = '~/Documents/MATLAB/COVIDProject/DataFilesUSCOVID-main/Plots'; 
%      saveas(gca, fullfile(fname, baseFileName), 'png');


TTT = sort(paramfitInts);
CILB = TTT(1+remove,:);
CIUB = TTT(end-remove,:);
mns = mean(paramfitInts);
meds = median(paramfitInts);
sdev = std(paramfitInts);
ARE = 100*mean(REs);
array2table(ARE)
array2table([CILB;CIUB])
r = [ARE;CILB;CIUB;mns;meds;sdev;paramfit];
end
end