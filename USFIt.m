function [] = USFIt
warning off;
close all
set(0,'DefaultFigureVisible','off')
% import data
%paramfits = readtable('paramfits.csv');
%paramfits = table2array(paramfits);
CasesByState = readtable('CasesByState.csv');
Dates = table2array(CasesByState(:,1));
CasesByState = table2array(CasesByState(:,2:end));
DeathsByState = readtable('DeathsByState.csv');
DeathsByState = table2array(DeathsByState(:,2:end));
TestingByState = readtable('CumTestsByState.csv');
TestingByState = table2array(TestingByState(:,2:end));
DatesByState = readtable('ImportantDates.csv');
DatesByState = table2array(DatesByState(:,2:end));
paramfits = readtable('params.csv');
paramfits = table2array(paramfits);

StateNames=["Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","District of Columbia","Florida","Georgia","Guam","Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Northern Mariana Islands","Ohio","Oklahoma","Oregon","Pennsylvania","Puerto Rico","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virgin Islands","Virginia","Washington","West Virginia","Wisconsin","Wyoming"];
Names = StateNames;
StatePops = [4903185,731545,7278717,3017804, 39512223, 5758736, 3565287, 973764, 705749,...
    21477737, 10617423,165768, 1415872, 1787065, 12671821, 6732219, 3155070, 2913314, 4467673,...
    4648794, 1344212, 6045680, 6892503, 9986857, 5639632, 2976149, 6137428, 1068778, 1934408,...
    3080156, 1359711, 8882190, 2096829, 19453561, 10488084, 762062,56882, 11689100, 3956971,...
    4217737,12801989,3193694,1059361, 5148714, 884659, 6829174, 28995881, 3205958, 623989,106977,...
    8535519, 7614893,1792147,5822434, 578759];
 T = 7.5;
 N = sum(StatePops);
 
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


for w = 1:258
    Cases(w) = sum(CasesByState(w,:));
    Deaths(w) = sum(DeathsByState(w,:));
    NewTests(w) = sum(TestingByState(w,:));
end
Cases = Cases';
Deaths = Deaths';

NewTests = NewTests';
NewTests = diff(NewTests);
    Stop = SHEs(52);
    alphBound = days(Stop-FirstCases(52))+21;
    FirstCase = FirstCases(52);
   

    Stop = days(Stop-FirstCase);
    SHE = [];
    TF = isnan(SHE);
    if TF == 1
        SHE = 59;
        SHS = 0;
    end
    SHE = Stop;
        
    lsqwt = Deaths(SHE+21)./Cases(SHE)

    
    
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
    
    %N = StatePops(its);
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
    k0=.005;

    % parameter bounds 
    Tg=6;
    betag=3/Tg;
    lsqwt=Deaths(SHE+21)/Cases(SHE);

   
  
   
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
   
        
    lb = 1;
    ub = 2000;
  
  %  paramguess = USparams(7)*USparams(1)*7.5;
    %fprintf('ITERATION %i OF %i\n',its,length(Names))
    tD = linspace(0,length(Cases)-1,length(Cases));
    
   
    %paramfit = [paramfit,21/paramfit(end)];
% Cfits = readtable('Cfits.csv');
% Cfits = table2array(Cfits);
% Dfits = readtable('Dfits.csv');
% Dfits = table2array(Dfits);
% URCfits = readtable('URCfits.csv');
% URCfits = table2array(URCfits);

%  for j = 1:Stop
%      Cfit(j) = sum(Cfits(j,:));
%      URCfit(j) = sum(URCfits(j,:));
%      TCT(j) = Cfit(j) + URCfit(j);
%  end
% 
%  for j = 1:Stop+21
%      Dfit(j) = sum(Dfits(j,:));
%  end


str = "United States";
%%%%%%%%%%%%%%%%%%%%%%%%
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
% tests = interp1(tD0,avs,tD);
    % [beta, psi, alpha, phi, t, tq, xi, k, A, I0, ag]
 %   lb = [0, 1, .005, .001,(max(tests(1:SHE))/(0.0050*N))*.8-eps,.4*max(tests(1:SHE)),1,.5,.005*N];
  %  ub = [6/7.5,6E3,1/7,.05,(max(tests(1:SHE))/(0.0050*N))*.8+eps,.6*max(tests(1:SHE)),800,12,.015*N];
% paramguess = (ub+lb)/2;
% [paramguess,resnorm] = lsqcurvefit(@objFxn2,paramguess,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
% array2table(paramguess)
% paramfit = paramguess;

% %        function [d1]=getData(CFit,nLevel,sets)
% %     d1 = zeros(sets,length(CFit));
% %     for k = 1:sets
% %         for j = 1:length(CFit)
% %             d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
% %             while d1(k,j) < 0
% %                 d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
% %             end
% %             d1(k,j) = ceil(d1(k,j));
% %         end
% %     end
% %    end
% %    nLevel = .5;
% %  
% %    sets = 10000;
% %    CdFit = diff(Cfit2);
% %    CdData = Cd(1:Stop-1);
% %    SynthC = getData(diff(Cfit2),nLevel,sets);
% %    figure
% % hold on
% % for jj = 1:10000
% % plot(1:1:130,SynthC(jj,:),'b','linewidth',2)
% % end
% % plot(1:1:130,CdFit,'k','linewidth',2)
% % plot(1:1:130,CdData,'r*','linewidth',2)
% % hold off
% % xlim([1,130])
% % SynthCumC = zeros(10000,131);
% % for jj = 1:10000
% %     for kk = 2:131
% %         if kk == 2
% %             SynthCumC(jj,kk-1) = Cases(1);
% %         end
% %         SynthCumC(jj,kk) = sum(SynthC(jj,1:kk-1));
% %     end
% % end
% % 
% % figure
% % hold on
% % for jj = 1:10000
% %     plot(0:1:130,SynthCumC(jj,:),'b','linewidth',2)
% % end
% % plot(0:1:130,Cfit2,'k','linewidth',2)
% % plot(0:1:130,Cases(1:Stop),'r*','linewidth',2)
% % hold off
% % xlim([0,130])
% % 
% % function z = objFxn2(x,tD)
% %         s0 = N;
% %         sq0 = x(9);
% %         k = x(5);
% %         A = x(6);
% %         z1 = zeros(length(tD),1);
% %         z2 = zeros(length(tD),1);
% %         %tests=movmean(NewTests,3);
% %         tests = interp1(tD0,avs,tD);
% %         rho = k*((tests/N))./(A/N+((tests/N)))+k0;
% %         mdd= zeros(length(tD),1);
% %         mid= zeros(length(tD),1);
% %         z1(1)=Cases(1);
% %         z2(1)=Deaths(1);
% % %         I0 = Cases(1)/rho(1);
% %         I0 = x(7);
% %         ag = x(8);
% %         bg = 21/ag;
% %         Iq0 = 0;
% %      
% %         %T = x(9);     
% %       
% %         u0 = [s0-sq0-I0,sq0,I0];
% %         xi = x(4);
% %         [t,Sol] = ode45(@(t,u)theModel(t,u,x),tD,u0);
% %         infs = Sol(:,3);
% %         infsq =  zeros(length(infs),1);
% %        % defs = Sol(:,5) + Sol(:,6);
% %        % z1= Sol(:,5); 
% %        % mass(1) = (gamcdf(1.5,5.1,.86)-gamcdf(0,5.1,.86))+(gamcdf(1.5,17.8,.45)-gamcdf(0,17.8,.45));
% %         for j = 2:length(tD)
% % %             mass(j) =(gamcdf(j+.5,5.1,.86)-gamcdf(j-.5,5.1,.86))+(gamcdf(j+.5,17.8,.45)-gamcdf(j-.5,17.8,.45));
% % %             defs(j) = mass(j)*defs(j);
% %             mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1);
% %             mdd(j-1)=(deathdist((2:j),ag,bg)-deathdist((1:j-1),ag,bg))*(infs(j-1:-1:1)+infsq(j-1:-1:1));
% %             z2(j) = z2(j-1)+mdd(j-1);
% %             z1(j) = z1(j-1)+rho(j-1)*mid(j-1);
% %             % get rid z2(j-1), z1(j-1)
% %         end
% % %         for j = 1:length(tD)
% % %             z1(j) = sum(infs(1:j));    
% % %         end
% %         z2 = z2*xi;
% %         z = [lsqwt*z1(1:SHE)',z2(1:SHE+21)'];
% % end
% % 
% % USparamfits = zeros(length(paramfit),10000);
% % USrhos = zeros(length(rho(1:Stop)),10000);
% % USR0 = zeros(length(rho(1:Stop)),10000);
% % USURCs = zeros(length(rho(1:Stop)),10000);
% % USC = zeros(length(rho(1:Stop)),10000);
% % USD = zeros(length(Deaths(1:Stop+21)),10000);
% % [paramguess,resnorm] = lsqcurvefit(@objFxn2,paramfit,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
% % [paramguess,resnorm] = lsqcurvefit(@objFxn2,paramguess,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
% % [paramguess,resnorm] = lsqcurvefit(@objFxn2,paramguess,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
% % save('paramfitUSA.mat','paramguess')
% % array2table(paramguess)
% % fprintf('---------------------------\n')
% % fprintf('STARTING CI GENERATION\N')
% % lb = [0, 10, .005, .001,.1905-eps,min(tests(75:SHE)),1,.5,.005*N];
% %     ub = [6/7.5,6E3,1/7,.05,.1905+eps,max(tests(75:SHE)),800,12,.05*N];
% % for jj = 1:10000
% %     [paramfit,resnorm] = lsqcurvefit(@objFxn2,paramguess,tD(1:SHE+21),[lsqwt*SynthCumC(jj,:),Deaths(1:SHE+21)'],lb,ub);
% %     USparamfits(:,jj) = paramfit;
% %     k = paramfit(5);
% %     A = paramfit(6);
% %     rho = k*((tests/N))./(A/N+((tests/N)))+k0;
% %     rho = rho(1:Stop);
% %     USrhos(:,jj) = rho;
% %     u0 = [N-paramfit(9)-paramfit(7),paramfit(9),paramfit(7)];
% %     [t,Sol] = ode45(@(t,u)theModel(t,u,paramfit),tD,u0);
% %     R0 = (Sol(:,1)/N)*7.5*paramfit(1);
% %     R0 = R0(1:Stop);
% %     USR0(:,jj) = R0;
% %     infs = Sol(:,3);
% %     infsq =  zeros(length(infs),1);
% %     z4(1) = paramfit(7)-Cases(1);
% %     for j = 2:length(tD(1:Stop))
% %         mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1);
% %     %mid(j-1)=infs(j-1)+infsq(j-1);
% %         z4(j) = z4(j-1)+(1-rho(j-1))*mid(j-1);
% %     end
% %     fits = objFxn2(paramfit,tD(1:SHE+21));
% %     USURCs(:,jj) = z4;
% %     USC(:,jj) = fits(1:SHE);
% %     USD(:,jj) = fits(SHE+1:end);
% %     if mod(jj,100) == 0
% %         fprintf('------------------------------\n')
% %         fprintf('GENERATING CIs %i PERCENT COMPLETE\n',round((jj/10000)*100))
% %         fprintf('------------------------------\n')
% %     end
% % end
% % save('USparamfits.mat','USparamfits');
% % save('USrhos.mat','USrhos');
% % save('USR0.mat','USR0') 
% % save('USURCs.mat','USURCs')
% % save('USC.mat','USC') 
% % save('USD.mat','USD')

   tests = interp1(tD0,avs,tD);
    % [beta, psi, alpha, phi, t, tq, xi, k, A, I0, ag]
    %lb = [0, 1, .005, .001,min(paramfits(:,5)),.4*max(tests(1:SHE)),1,.5,.0025*N];
    %ub = [6/7.5,1000,1/7,.05,.525*max(paramfits(:,5)),.6*max(tests(1:SHE)),800,12,.015*N];
   
 
    %if its == 12
    %    lb(6) = .15*max(NewTests(1:SHE+21));
    %    ub(6) = .35*max(NewTests(1:SHE+21));
    %end
    paramguess = [0.286338928636279,258.887865817671,0.0240070128944282,0.00861721915553062,0.358215252628226,271068.991836735,50.3949836983015,12.0000000000000,829407.113710649,1.75000000000000,3.50535973190290,0.0367698558305923,0.0990264973789266,0.0101613689219439];
    paramguess = paramguess(1:9);
    lb = paramguess - eps*ones(1,9);
    ub = paramguess + eps*ones(1,9);
    %if its == 53
    %    paramguess = .5*(ub-lb);
    %end
    [paramfit,resnorm] = lsqcurvefit(@objFxn2,paramguess,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
    [paramfit,resnorm] = lsqcurvefit(@objFxn2,paramfit,tD(1:SHE+21),[lsqwt*Cases(1:SHE)',Deaths(1:SHE+21)'],lb,ub);
    paramfit = [paramfit,21/paramfit(end-1)];


    fits = objFxn2(paramfit,tD);
    Cfit = fits(1:SHE);
    Dfit = fits(SHE+1:end);
    save('fits2.mat','fits')
    k = paramfit(5);
    A = paramfit(6);
    rho = k*((tests/N))./(A/N+((tests/N)))+k0;
     


%str = Names(its);
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
%title(str)
legend('Cumulative Death Data','Model Fit','Stay at Home Order in Effect','FontSize',10,'Location','NorthWest')
hold off
set(gca,'FontSize',16)
%baseFileName = sprintf('CumDeaths%d',its);
%fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
%saveas(gca, fullfile(fname, baseFileName), 'png');
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
figure
plot(tD(1:SHE),rho(1:SHE),'b','linewidth',2)
ylabel('Capture Rate')
xlabel('Days since First Reported Case up to May 31')
ylim([0,1.05*max(rho(1:SHE))])
xlim([tD(1), tD(SHE)])
set(gca,'FontSize',16)

%   baseFileName = sprintf('rho%d',its);
%    fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
%    saveas(gca, fullfile(fname, baseFileName), 'png');

pos = (movmean(Cd,3)./tests(1:end-1))*100;
posFit =(diff(Cfit/lsqwt)'./tests(1:SHE-1))*100;
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
   legend({'Positives Inferred From Data','Model Predicted Positives','CDC Positivity Threshold','New Tests Interpolaton','New Daily Tests Raw Data'},'FontSize',10,'Location','North')
   ylabel('New Daily Tests (% population)')
   xlabel('Days Since First Reported Case up to May 31')
   xlim([tP(1) tP(end)])
   xticks(round(linspace(tP(1),tP(end),9)))
   yticks(round(linspace(0,1.05*max(((tests(tP+2)./N)*100)),6),2))
   ylim([0,1.05*max(((tests(tP+2)./N)*100))])
   %title(str)
   set(gca,'FontSize',16)
 %   baseFileName = sprintf('Positives%d',its);
 %   fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
 %   saveas(gca, fullfile(fname, baseFileName), 'png');

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

u0 = [N-paramguess(9)-paramguess(7),paramguess(9),paramguess(7)];
[t,Sol] = ode45(@(t,u)theModel(t,u,paramfit),tD,u0);
infs = Sol(:,3);
R0 = (Sol(:,1)/N)*7.5*paramguess(1);
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
%title(str)
legend('Reported Cumulative Case Data','Fit','Unreported Case Estimate','True Cumulative Case Estimate','Stay at Home Order in Effect','FontSize',10,'Location','NorthWest')
hold off
set(gca,'FontSize',16)
%baseFileName = sprintf('CumCases%d',its);
%fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
%saveas(gca, fullfile(fname, baseFileName), 'png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TCT = (z3(1:SHE))+((Cfit))/lsqwt;
fprintf('Case Ratio\n')
TCT(end)/(Cfit(end)/lsqwt)
%TCTd = diff(TCT);

%figuregca
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
figure
hold on
plot(tD(2:SHE),Cd(1:SHE-1),'r*',tD(2:SHE),diff(Cfit)/lsqwt,'b','LineWidth',2)
plot(tD(2:SHE),diff(TCT),'b-.','LineWidth',2)
hold off
ylim([0,1.05*max(diff(TCT))])
xlim([tD(2),tD(SHE)])
xlabel('Days Since First Reported Case up to May 31')
ylabel('Count')
set(gca,'FontSize',16)
paramfit = [paramfit,(TCT(end)/N)*100,(Dfit(end)/N)*100,mean(rho(1:SHE)),Cfit(end)/lsqwt/TCT(end)];
save('USAparamCurr.mat','paramfit')
figure
hold on
plot(tD(1:SHE),R0(1:SHE),'b','LineWidth',2)
plot(tD(1:SHE),ones(1,SHE),'k-.')
hold off
ylabel('Reproduction Number')
xlabel('Days Since First Reported Case up to May 31')
ylim([.95*min(R0(1:SHE)),1.05*max(R0(1:SHE))])
xlim([tD(1),tD(SHE)])
set(gca,'FontSize',16)
array2table(paramfit)
tests = tests(1:SHE);
URC = z3(1:SHE);
save('URCUSACurr.mat','URC')
save('testsUSACurr.mat','tests')
save('CfitUSACurr.mat','Cfit')
save('TCTUSACurr.mat','TCT')
rho = rho(1:SHE);
save('rhoUSACurr.mat','rho')
save('R0USACurr.mat','R0')

paramguess = paramfit(1:9);
lb = .25*paramguess;
     ub = 4*paramguess;
     lb(5) = .8*paramguess(5);
     ub(5) = 1.2*paramguess(5);
     lb(1) = 0;
     ub(1) = 6/7.5;
     lb(7) = max(lb(7),1);
     lb(6) = .4*max(tests(1:SHE));
     ub(6) = .6*max(tests(1:SHE));
     lb(9) = .0025*N;
     ub(9) = .015*N;
paramguess2 = (ub+lb)./2;



       function [d1]=getData(CFit,nLevel,sets)
    d1 = zeros(sets,length(CFit));
    for k = 1:sets
        for j = 1:length(CFit)
            d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
            while d1(k,j) < 0
                d1(k,j) = CFit(j) + normrnd(0,CFit(j)*nLevel.^2);
            end
            d1(k,j) = ceil(d1(k,j));
        end
    end
   end
   nLevel = .6;
 
   sets = 10000;
   CdFit = diff(Cfit)/lsqwt;
   CdData = Cd(1:Stop-1);
   SynthC = getData(diff(Cfit)/lsqwt,nLevel,sets);
   fprintf('Generating plot 1\n')
   figure
hold on
for jj = 1:sets
plot(1:1:130,SynthC(jj,:),'b','linewidth',2)
end
plot(1:1:130,CdData,'r*','linewidth',1)
plot(1:1:130,CdFit,'k','linewidth',2)
hold off
xlim([1,130])
SynthCumC = zeros(sets,131);
for jj = 1:sets
    for kk = 2:131
        if kk == 2
            SynthCumC(jj,kk-1) = Cases(1);
        end
        SynthCumC(jj,kk) = sum(SynthC(jj,1:kk-1));
    end
end
baseFileName = sprintf('SyntheticDailyCases');
fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');
fprintf('generating plot 2')
figure
hold on
for jj = 1:sets
    plot(0:1:130,SynthCumC,'b','linewidth',2)
end
plot(0:1:130,Cases(1:131),'r*','linewidth',1)
plot(0:1:130,Cfit/lsqwt,'k','linewidth',2)
hold off
xlim([0,130])
baseFileName = sprintf('SyntheticCumCases');
fname = 'C:\Users\macdo\OneDrive\Desktop\TestingFit\Plots';
saveas(gca, fullfile(fname, baseFileName), 'png');
% 
% USparamfits = zeros(9,sets);
% USrhos = zeros(length(rho(1:Stop)),sets);
% USR0 = zeros(length(rho(1:Stop)),sets);
% USURCs = zeros(length(rho(1:Stop)),sets);
% USC = zeros(length(rho(1:Stop)),sets);
% USD = zeros(length(Deaths(1:Stop+21)),sets);
% 
% fprintf('STARTING CI GENERATION\N')
% %paramguess = paramfit(1:9);
% clear paramfit
% for jj = 1:sets
%     lsqwt = Deaths(SHE+21)/SynthCumC(jj,end);
%     lb = .25*paramguess;
%     ub = 4*paramguess;
%     lb(5) = .8*paramguess(5);
%     ub(5) = 1.2*paramguess(5);
%     lb(1) = 0;
%     ub(1) = 6/7.5;
%     lb(7) = max(lb(7),1);
%     lb(6) = .4*max(tests(1:SHE));
%     ub(6) = .6*max(tests(1:SHE));
%     lb(9) = .0025*N;
%     ub(9) = .015*N;
%     [paramfit,resnorm] = lsqcurvefit(@objFxn2,paramguess,tD(1:SHE+21),[lsqwt*SynthCumC(jj,:),Deaths(1:SHE+21)'],lb,ub);
%     USparamfits(:,jj) = paramfit;
%     array2table(paramfit)
%     k = paramfit(5);
%     A = paramfit(6);
%     rho = k*((tests/N))./(A/N+((tests/N)))+k0;
%     rho = rho(1:Stop);
%     USrhos(:,jj) = rho;
%     u0 = [N-paramfit(9)-paramfit(7),paramfit(9),paramfit(7)];
%     [t,Sol] = ode45(@(t,u)theModel(t,u,paramfit),tD,u0);
%     R0 = (Sol(:,1)/N)*7.5*paramfit(1);
%     R0 = R0(1:Stop);
%     USR0(:,jj) = R0;
%     infs = Sol(:,3);
%     infsq =  zeros(length(infs),1);
%     z4(1) = paramfit(7)-Cases(1);
%     for j = 2:length(tD(1:Stop))
%         mid(j-1)=(infdist((2:j),1/T)-infdist((1:j-1),1/T))*infs(j-1:-1:1);
%     %mid(j-1)=infs(j-1)+infsq(j-1);
%         z4(j) = z4(j-1)+(1-rho(j-1))*mid(j-1);
%     end
%     fits = objFxn2(paramfit,tD(1:SHE+21));
%     USURCs(:,jj) = z4;
%     USC(:,jj) = fits(1:SHE)/lsqwt;
%     USD(:,jj) = fits(SHE+1:end);
%     if mod(jj,100) == 0
%         fprintf('------------------------------\n')
%         fprintf('GENERATING CIs %i PERCENT COMPLETE\n',round((jj/sets)*100))
%         fprintf('------------------------------\n')
%     end
% end
% save('USparamfits.mat','USparamfits');
% save('USrhos.mat','USrhos');
% save('USR0.mat','USR0') 
% save('USURCs.mat','USURCs')
% save('USC.mat','USC') 
% save('USD.mat','USD')

end