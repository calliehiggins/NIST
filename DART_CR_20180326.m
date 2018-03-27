%Determining the resonant frequency, Q, Adrive, and phidrive using DART
%atomic force acoustic microscopy, model used from A. Gannepalli et al
%Nanotechnology 22 (2011) paper (with corrigendum)
%Created by Callie Higgins, 15 November 2017
clear all
%-------------------------------------------------------------------------
% Setting scaling factors
%-------------------------------------------------------------------------
giga = 10^9;
mega = 10^6;
kilo = 10^3;
deci = 10^(-1);
centi = 10^(-2); 
milli = 10^(-3);
micro = 10^(-6);
nano = 10^(-9);

%Importing relevant variables (phase, amplitude, and frequency for both
%lockin amplifiers

%Importing all Dual Cure Files from all cubes

cube01= csvread('VariedPowerAtTipCRv2.csv',2,0); %not really Cube01, just didn't want to go through whole code and change
crTune= csvread('CRTune.csv',2,0);
%cube01= csvread('Power4mW_4distances.csv',2,0); %not really Cube01, just didn't want to go through whole code and change
%cube01= csvread('Power8mW_4distances.csv',2,0); %not really Cube01, just didn't want to go through whole code and change

%cube01= csvread('dC01Scans.csv',2,0);
%cube10= csvread('fC10Scans.csv',2,0);
%cube30= csvread('gC30Scans.csv',2,0);
% %resin= csvread('ResinCRTuneC1LasD0003.csv',2,0);
% %resin= csvread('ResinCRTuneC1LasE0004.csv',2,0);
%resin= csvread('ResinCRTuneC1LasE0006.csv',2,0);
%resin= csvread('CR_AuGlyFC0002.csv',2,0); %Actually Glycerol

%% Assigned input variables for all (adapted from kvf_lateral_viscoelastic.nb code)
QFR=100;%free res Q
fFR=13.2*kilo; %free res resonant freq
% xL0= 1.8751041; %for mode one
% tang = 12*pi./180; % Tip angle
% hL = 0.05; %tip height h divided by length of cantilever to tip (could be 0.02)
% %  r = [0.85 ,0.9,0.95,1];    %ratio klat/kvert (how much lateral dashpot/spring is contributing to measurement
% %  gamma = [0.9,0.95, 1, 1.05]; %L1/L0 (location of tip wrt length of cantilever)
% r = 1;    %ratio klat/kvert (how much lateral dashpot/spring is contributing to measurement
gamma = 1.035; %L1/L0 (location of tip wrt length of cantilever)
L=240*nano; %Length of Cantilever in meters
kc=3.3173; %Spring constant of cantilever N/m

% sa= sin(tang);
% ca= cos(tang);
% 
% c2rs2 = ca^2 + r.*sa^2;
% s2rc2 = sa^2 + r.*ca^2;
% cs1r = sa.*ca.*(r - 1);
% pre = 0.5*r.*hL.^2;
% 
% pos = [0,1,5,10,50]; %positions of the laser wrt cantilever tip

% Extracting f1, f2, A1, A2, phi1, and phi2 of all curves 
[t01, A101, A201, phi101, phi201, f101, f201] = DARTCRVars(cube01,5,30000,1);
[f0n, A0n] =TuneVars(crTune,1,1);
fitParm0=2200;
fitParm1=4000;
fn0(:,:)=f0n{1}(fitParm0:fitParm1,1);
A0n0(:,:)=A0n{1}(fitParm0:fitParm1,1);
%dPhi= 10:5:30;
dPhi=25;
b1=150000;
avg=100;
f=(b1:500:b1+300000);

i0=[1,2,4,8,12];
%i0=[4,4,4,4];
%i0=[8,8,8,8];

% d0{1}=i0(1).*(t01{1});%-1.5);
% d0{2}=i0(2).*(t01{2});%-1.5);%subtraction 1.5 for varied I, 5.5 for varied distance stat (2-5)
% d0{3}=i0(3).*(t01{3});%-1.5);
% d0{4}=i0(4).*(t01{4});%-1.5);
% d0{5}=i0(5).*(t01{5});%1.5);
% %d0{5}=i0(5).*(t01{5}-1.5): %for varied exposure intensity

d0{1}=i0(1).*(t01{1}-0.044);
d0{2}=i0(2).*(t01{2}-0.065);%subtraction 1.5 for varied I, 5.5 for varied distance stat (2-5)
d0{3}=i0(3).*(t01{3}-0.065);
d0{4}=i0(4).*(t01{4}-0.065);
d0{5}=i0(5).*(t01{5}-0.26);


%% Variable Power at tip
%Plotting different delta phi
% [f0_0,Q_0,Ad_0,phid_0,omega0,phi0,X10,X20,An0,fmean0,Admean0,Qmean0]= solvCR(t01{1}, A101{1}, A201{1}, phi101{1}, phi201{1}, f101{1}, f201{1},dPhi(1,1), avg,f);
% [f0_1,Q_1,Ad_1,phid_1,omega1,phi1,X11,X21,An1,fmean1,Admean1,Qmean1]= solvCR(t01{1}, A101{1}, A201{1}, phi101{1}, phi201{1}, f101{1}, f201{1},dPhi(1,2), avg,f);
% [f0_2,Q_2,Ad_2,phid_2,omega2,phi2,X12,X22,An2,fmean2,Admean2,Qmean2]= solvCR(t01{1}, A101{1}, A201{1}, phi101{1}, phi201{1}, f101{1}, f201{1},dPhi(1,3), avg,f);
% [f0_3,Q_3,Ad_3,phid_3,omega3,phi3,X13,X23,An3,fmean3,Admean3,Qmean3]= solvCR(t01{1}, A101{1}, A201{1}, phi101{1}, phi201{1}, f101{1}, f201{1},dPhi(1,4), avg,f);
% [f0_4,Q_4,Ad_4,phid_4,omega4,phi4,X14,X24,An4,fmean4,Admean4,Qmean4]= solvCR(t01{1}, A101{1}, A201{1}, phi101{1}, phi201{1}, f101{1}, f201{1},dPhi(1,5), avg,f);
% plot(f./1000,An0)
% hold on
% plot(f./1000,An1)
% plot(f./1000,An2)
% plot(f./1000,An3)
% plot(f./1000,An4)
% plot(fn0./1000,A0n0)
% xlabel('Frequency (kHz)','fontsize',15)
% ylabel('Amplitude (V)','fontsize',15)
% legend(gca,'\Delta\phi = 10','\Delta\phi = 15','\Delta\phi = 20',...
%     '\Delta\phi = 25','\Delta\phi = 30','CR Tune')
% %legend(gca,'0 \mum','1 \mum','5 \mum','10 \mum')
% set(gcf,'color','w');
% set(gca, 'YScale', 'log')
% %xlim([10^(-2) 450])
% box on
% grid on
% hold off

%% Calculating SHO f0, Q, Ad, and Phid
[f0_01{1},Q_01{1},Ad_01{1},phid_01{1},omega{1},phi{1},X1{1},X2{1},An{1},fmean{1},Admean{1},Qmean{1}]= solvCR(t01{1}, A101{1}, A201{1}, phi101{1}, phi201{1}, f101{1}, f201{1},dPhi(1,1),avg,f);
[f0_01{2},Q_01{2},Ad_01{2},phid_01{2},omega{2},phi{2},X1{2},X2{2},An{2},fmean{2},Admean{2},Qmean{2}]= solvCR(t01{2}, A101{2}, A201{2}, phi101{2}, phi201{2}, f101{2}, f201{2},dPhi(1,1),avg,f);
[f0_01{3},Q_01{3},Ad_01{3},phid_01{3},omega{3},phi{3},X1{3},X2{3},An{3},fmean{3},Admean{3},Qmean{3}]= solvCR(t01{3}, A101{3}, A201{3}, phi101{3}, phi201{3}, f101{3}, f201{3},dPhi(1,1),avg,f);
[f0_01{4},Q_01{4},Ad_01{4},phid_01{4},omega{4},phi{4},X1{4},X2{4},An{4},fmean{4},Admean{4},Qmean{4}]= solvCR(t01{4}, A101{4}, A201{4}, phi101{4}, phi201{4}, f101{4}, f201{4},dPhi(1,1),avg,f);
[f0_01{5},Q_01{5},Ad_01{5},phid_01{5},omega{5},phi{5},X1{5},X2{5},An{5},fmean{5},Admean{5},Qmean{5}]= solvCR(t01{5}, A101{5}, A201{5}, phi101{5}, phi201{5}, f101{5}, f201{5},dPhi(1,1),avg,f);
%plot(dPhi, Q_01{2}(1,:))
% xlabel('\Delta \phi','fontsize',15)
% ylabel('Q','fontsize',15)
%legend(gca,'Lorentz Fit of DART @ t=0','CR Tune')
%legend(gca,'0 \mum','1 \mum','5 \mum','10 \mum')
% set(gcf,'color','w');
% %set(gca, 'YScale', 'log')
% %xlim([10^(-2) 450])
% box on
% grid on
% hold off
%% Plotting CR Tune Raw Data vs. Lorentz Fit of DART data
% b1=150000;
% avg=1;
% f=(b1:500:b1+300000);
% f0_n=mean(f0_01{1}(1:avg,1));
% Ad_n=mean(Ad_01{1}(1:avg,1));
% Q_n=mean(Q_01{1}(1:avg,1));
% An=f0_n.^2.*Ad_n./((f0_n.^2-f.^2).^2+(f0_n.*f./Q_n).^2).^(1/2);
% plot(f./1000,An)
% %hold on
% %plot(f0n{1}./1000,A0n{1})
% %plot(f0_01{1}./1000,A101{1})
% %plot(f0_01{1}./1000,A201{1})
% xlabel('Frequency (kHz)','fontsize',15)
% ylabel('Amplitude (V)','fontsize',15)
% legend(gca,'Lorentz Fit of DART @ t=0','CR Tune')
% %legend(gca,'0 \mum','1 \mum','5 \mum','10 \mum')
% set(gcf,'color','w');
% set(gca, 'YScale', 'log')
% %xlim([10^(-2) 450])
% box on
% grid on
% hold off

%% Calculating the loss tangent

%Without Lateral spring and dashpot
[chi_01,lL_01,s_01,alpha_01,beta_01,ks_01, lossT_01,z_01,a_01,b_01,alphaR_01,ksR_01,lossTR_01]=ElasticCalc(QFR,fFR,Q_01,f0_01,kc,gamma,L);

% With Lateral spring and dashpot can only be calculated with singe value for gamma and r
%[chi_01L,CharN_01L,CharP_01L, alpha_01L,beta_01L,ks_01L, lossT_01L,alphaR_01L,betaR_01L,ksR_01L,lossTR_01L,z_01L]=ECalcLatSND(QFR,fFR,Q_01,f0_01,xL0,gamma,c2rs2,s2rc2,cs1r,pre,hL);


%% Plotting Data
%figure(6)
%plotFQPA(t01,f0_01,Q_01, phid_01, Ad_01, A101, A201, phi101, phi201, f101, f201)

x=t01; %d0 or t01
xTitle = 'Time (s)';%'Energy Dose (mJ)' or 'Time (s)'
plotVar(x,Q_01,xTitle,'Q','',1)
plotVar(x,f0_01,xTitle,'f_0 (kHz)','',2)
plotVar(x,phid_01,xTitle,'\phi_{drive}','',3)
plotVar(x,Ad_01,xTitle,'Amplitude_{drive} (V)','',4)
plotVar(x,lossT_01,xTitle,'tan\delta_{tip-sample}','',5)
%plotVar(x,lossTR_01,xTitle,'tan\delta_{tip-sample_f} / tan\delta_{tip-sample_i} ','',6)
plotVar(x,alpha_01,xTitle,'\alpha, Contact Stiffness','',7)
%plotVar(x,alphaR_01,xTitle,'\alpha_f / \alpha_i','',8)
plotVar(x,ks_01,xTitle,'k_{system} (N/m)','',9)
%plotVar(x,ksR_01,xTitle,'k_{system_f} / k_{system_i}','',10)
% plotVar(t01,ks_01,xTitle','k_s_y_s_t_e_m','',8)
% plotVar(d0,beta_01,'Energy Dose (mJ))','\beta_{tip-sample}, Damping Coefficient','',8)
% plotVar(t01,a_01,xTitle,'a','',9)
% plotVar(t01,b_01,xTitle,'b','',10)
% plotVar(t01,b_01,xTitle,'s','',11)




%% Function definitions %%
function [t, A1, A2, phi1, phi2, f1, f2]=DARTCRVars(cubeName, varNum, deltaf,tVar)    %varNum indicates number of different location tested
    for i=1:1:size(cubeName, 2)
        for j = 1:1:varNum
                c = (j-1)*6+1;
                t{j}=cubeName(tVar:end-10,c)-1;
                d = (j-1)*6+2;
                A1{j}=(cubeName(tVar:end-10,d));
                e = (j-1)*6+3;
                A2{j}=cubeName(tVar:end-10,e);
                f = (j-1)*6+4;
                phi1{j}=cubeName(tVar:end-10,f);
                g = (j-1)*6+5;
                phi2{j}=cubeName(tVar:end-10,g);
                h = (j-1)*6+6;
                f1{j}=cubeName(tVar:end-10,h); %for when f=f1
                f2{j}=f1{j}+deltaf;
%                 f1{j}=cubeName(:,h)-deltaf/2; %for when the f=fc
%                 f2{j}=f1{j}+deltaf; 
%                 f2{j}=cubeName(:,h);
%                 f1{j}=f2{j}-deltaf; %for when f=f2
        end
    end
end

function [f, A1]=TuneVars(cubeName, varNum, tVar)    %varNum indicates number of different location tested
    for i=1:1:size(cubeName, 2)
        for j = 1:1:varNum
                c = (j-1)*2+1;
                f{j}=cubeName(tVar:end-10,c)-1;
                d = (j-1)*2+2;
                A1{j}=(cubeName(tVar:end-10,d));
        end
    end
end


%[f0,Q,Ad,phid,omega,phi,X1,X2,f, An,f0_n,Q_n]
function [f0,Q,Ad,phid,omega,phi,X1,X2,An,fmean,Admean,Qmean]= solvCR(t, A10, A20, phi1, phi2, f1, f2,deltaPhi, avg,f)
        A1= A10;%-(-3.41757*10^(-8)*f1+0.0049);
        A2= A20;%-(-3.41757*10^(-8)*f2+0.0049);
        omega = (f1.*A1)./(f2.*A2);
        %phi = tan(degtorad(phi2-phi1-25));
        phi = tan(degtorad(phi2-phi1+deltaPhi));
        X1= -(1-(sign(phi).*(1+phi.^2).^(1/2)./omega))./phi;
        X2= (1-sign(phi).*omega.*(1+phi.^2).^(1/2))./phi;
        %f0 =(sqrt(f1.*f2.*(f2.*X1-f1.*X2)./(f1.*X1-f2.*X2)));
        f0 =sqrt(f1.*f2.*((f2.*X1-f1.*X2)./(f1.*X1-f2.*X2)));
        %f0 = (f0-mean(f0(1:250))+246900);

        Q = sqrt((f1.*f2.*(f2.*X1-f1.*X2).*(f1.*X1-f2.*X2)))./(f2.^2-f1.^2);
        %Q = Q0-mean(Q0(100:1000,1))+7;
        Ad=A1.*sqrt((f0.^2-f1.^2).^2+(f0.*f1./Q).^2)./f0.^2;
        %Ad=Ad0-mean(Ad0(100:1000,1))+1.6*10^(-3);

        phid = phi1 - radtodeg(atan((f0.*f1)./(Q.*(f0.^2-f1.^2))));%+12.5);
        fmean=mean(f0(1:avg,1));
        Admean=mean(Ad(1:avg,1));
        Qmean=mean(Q(1:avg,1));
        An=fmean.^2.*Admean./((fmean.^2-f.^2).^2+(fmean.*f./Qmean).^2).^(1/2);
end

% New CR solver with variable delta phi incorporated
function [f0,Q,Ad,phid,omega,phi,X1,X2]= solvCRn(t, A10, A20, phi1, phi2, f1, f2,deltaPhi)
        A1= A10;%-(-3.41757*10^(-8)*f1+0.0049);
        A2= A20;%-(-3.41757*10^(-8)*f2+0.0049);
        omega = (f1.*A1)./(f2.*A2);
        phi = tan(degtorad(phi2-phi1+deltaPhi));
        X1= -(1-(sign(phi).*(1+phi.^2).^(1/2)./omega))./phi;
        X2= (1-sign(phi).*omega.*(1+phi.^2).^(1/2))./phi;
        f0 =sqrt(f1.*f2.*((f2.*X1-f1.*X2)./(f1.*X1-f2.*X2)));
        Q = sqrt((f1.*f2.*(f2.*X1-f1.*X2).*(f1.*X1-f2.*X2)))./(f2.^2-f1.^2);
        Ad=A1.*sqrt((f0.^2-f1.^2).^2+(f0.*f1./Q))./f0.^2;
        phid = phi1 - radtodeg(atan((f0.*f1)./(Q.*(f0.^2-f1.^2))));%+12.5);
end

function plotFQPA(t, f0, Q, phid, Ad, A1, A2, phi1, phi2, f1, f2)  
    sz = 1;
    for i=1:length(t) 
        %u=length(t);
        figure(i)
        subplot(4,1,1)
        plot(t{i},f1{i})
        hold on
        plot(t{i},f2{i})
        scatter(t{i},f0{i},sz,'.k')
        ylabel('f (Hz)')
        legend(gca,'f_1','f_2','f_0')
        %hold off

        subplot(4,1,2)
        scatter(t{i},Q{i},sz,'.k')
        ylabel('Q')

        subplot(4,1,3)
        plot(t{i},phi1{i})
        hold on
        plot(t{i},phi2{i})
        scatter(t{i},phid{i},sz,'.k')
        hold off
        ylabel('\phi')
        legend(gca,'\phi_1','\phi_2','\phi_D')

        subplot(4,1,4)
        plot(t{i},A1{i})
        hold on
        plot(t{i},A2{i})
        scatter(t{i},Ad{i},sz,'.k')
        ylabel('Amplitude (V)')
        xlabel('Time (s)')
        legend(gca,'A_1','A_2','A_D')
    end
end

function [chi,lL,s,alpha,beta,ks,lossT, z,a,b,alphaR,ksR,lossTR]=ElasticCalc(QFR,fFR,Q,f0,kc,gamma,L)
%from Killgore Langmuir 2011
    for j= 1
        beg0=1900;
        en0=2200;
        %Q=100;
        chi{j} = 2*pi*fFR./QFR;
        xL = 1.8751; %for the first free flexural vibration mode
        a{j} = xL.*(f0{j}./fFR).^(1/2);
        b{j} = a{j}.*((2*pi.*f0{j}-chi{j}.*Q{j})./(8*pi.*f0{j}.*Q{j}));
        lL{j} = a{j}+i.*b{j};
        %lL1{j} = lL{j}.*gamma; %because the tip is at the end of the cantilever, L=L1
        %Lp = L.*(1-gamma); %Lp=L-L1, L=length of cantiliver, L1 = location of tip on cant
        s{j} = (2/3).*(lL{j}.*gamma).^3.*(1+cos(lL{j}).*cosh(lL{j}))./...
            ((sinh(lL{j}.*gamma).*cos(lL{j}.*gamma)-sin(lL{j}.*gamma).*cosh(lL{j}.*gamma)).*...
            (1+cos(lL{j}.*(1-gamma)).*cosh(lL{j}.*(1-gamma)))+...
            (1-cos(lL{j}.*gamma).*cosh(lL{j}.*gamma)).*...
            (sin(lL{j}.*(1-gamma)).*cosh(lL{j}.*(1-gamma))-cos(lL{j}.*(1-gamma)).*sinh(lL{j}.*(1-gamma))));
        alpha{j}=real(s{j});
        %alpha{j}=alpha0{j}-mean(alpha0{j}(1900:1950));
        alphaR{j}=alpha{j}./(mean(alpha{j}(beg0:en0)));
        ks{j}=alpha{j}.*kc;
        ksR{j}=ks{j}./(mean(ks{j}(beg0:en0)));
        beta{j}=imag(s{j})./((a{j}.^2-b{j}.^2));
        lossT{j} = (beta{j}.*f0{j}.*xL.^2)./(alpha{j}.*fFR);
        lossTR{j}=lossT{j}./(mean(lossT{j}(beg0:en0)));
        z{j} = 1:1:length(lossT{j});
    end
    for j=2:length(Q)
        beg=0;%1900;
        en=0;%2200;
        %Q=100;0
        chi{j} = 2*pi*fFR./QFR;
        xL = 1.8751; %for the first free flexural vibration mode
        a{j} = xL.*(f0{j}./fFR).^(1/2);
        b{j} = a{j}.*((2*pi.*f0{j}-chi{j}.*Q{j})./(8*pi.*f0{j}.*Q{j}));
        lL{j} = a{j}+i.*b{j};
        %lL1{j} = lL{j}.*gamma; %because the tip is at the end of the cantilever, L=L1
        %Lp = L.*(1-gamma); %Lp=L-L1, L=length of cantiliver, L1 = location of tip on cant
        s{j} = (2/3).*(lL{j}.*gamma).^3.*(1+cos(lL{j}).*cosh(lL{j}))./...
            ((sinh(lL{j}.*gamma).*cos(lL{j}.*gamma)-sin(lL{j}.*gamma).*cosh(lL{j}.*gamma)).*...
            (1+cos(lL{j}.*(1-gamma)).*cosh(lL{j}.*(1-gamma)))+...
            (1-cos(lL{j}.*gamma).*cosh(lL{j}.*gamma)).*...
            (sin(lL{j}.*(1-gamma)).*cosh(lL{j}.*(1-gamma))-cos(lL{j}.*(1-gamma)).*sinh(lL{j}.*(1-gamma))));
        alpha{j}=real(s{j});
        %alpha{j}=alpha0{j}-mean(alpha0{j}(5700:5800));
        %alphaR{j}=alpha{j}./(mean(alpha{j}(beg:en)));
        ks{j}=alpha{j}.*kc;
        %ksR{j}=ks{j}./(mean(ks{j}(beg:en)));
        beta{j}=imag(s{j})./((a{j}.^2-b{j}.^2));
        lossT{j} = (beta{j}.*f0{j}.*xL.^2)./(alpha{j}.*fFR);
        %lossTR{j}=lossT{j}./(mean(lossT{j}(beg:en)));
        z{j} = 1:1:length(lossT{j});
    end
end
% 
% function [B]=curvefit(Var1, Var2)
%    for j=1:length(Var1)     
%         F{j} = @(B,Var1) min(B(1),B(2)+B(3)*Var1);   %Form of the equation
%         IC = [max(Var2{j}) max(Var2{j}) 0]; %Initial guess
%         B{j} = lsqcurvefit(F,IC,Var1{j},Var2{j},[min(Var2{j}) -inf -inf],[max(Var1{j}) inf 0]);
%         a = (B{j}(1) - B{j}(2)) / B{j}(3);
%         cte = B{j}(1);
%         c = B{j}(2);
%         d = B{j}(3);
%    end
% end

function [chi,CharN,CharP, alpha,beta,ks,lossT,alphaR,betaR,ksR,lossTR,z]=ECalcLatSND(QFR,fFR,Q,f0,xL0,gamma,c2rs2,s2rc2,cs1r,pre,hL) %elastic calculations incorporating the lateral spring and dashpot compontents
    for j=1:length(Q)  
        chi = 2*pi*fFR./QFR;
        a = xL0.*(f0{j}./fFR).^(1/2);
        b = a.*((2*pi.*f0{j}-chi.*Q{j})./(8*pi.*f0{j}.*Q{j}));
        xLna = a + b.*1i;
        xLn1a = gamma.*xLna;
        xLnpa = (1-gamma).*(xLna);
        Aa = pre.*(1 - cos(xLn1a).*cosh(xLn1a)).*(1 + cos(xLnpa).*cosh(xLnpa));
        Ca = (xLn1a.^4).*(1 + cos(xLna).*cosh(xLna));
        B1a = 0.5*hL.^2.*xLn1a.^3.*s2rc2.*(((1 + cos(xLnpa).*cosh(xLnpa)).*(sin(xLn1a) ...
            .*cosh(xLn1a) + cos(xLn1a).*sinh(xLn1a))) - (1 -cos(xLn1a).*cosh(xLn1a))...
            .*(sin(xLnpa).*cosh(xLnpa) + cos(xLnpa).*sinh(xLnpa)));
        B2a = hL.*xLn1a.^2.*(cs1r.*(((1 + cos(xLnpa).*cosh(xLnpa)).*sin(xLn1a).*...
            sinh(xLn1a)) + (1 - cos(xLn1a).*cosh(xLn1a)).*sin(xLnpa).*sinh(xLnpa)));
        B3a = 0.5.*xLn1a.*c2rs2.*(((1 + cos(xLnpa).*cosh(xLnpa)).*(sin(xLn1a).*...
            cosh(xLn1a) - cos(xLn1a).*sinh(xLn1a)) - (1 - cos(xLn1a).*cosh(xLn1a)).*...
           (sin(xLnpa).*cosh(xLnpa) - cos(xLnpa).*sinh(xLnpa))));
        Ba = B1a+B2a+B3a;
        CharP{j} = (-Ba + sqrt(Ba.^2 - 4.*Aa.*Ca))./(6 .*Aa);
        CharN{j} = (-Ba - sqrt(Ba.^2 - 4.*Aa.*Ca))./(6 .*Aa);
        beta{j} = imag(CharP{j}./((a.^2-b.^2)).*gamma.^2);
        betaR(1,j) = mean(beta{j}(4500:end,1))./mean(beta{j}(100:500,1));        
        alpha{j} = real(CharP{j});
        alphaR(1,j) = mean(alpha{j}(4500:end,1))./mean(alpha{j}(100:500,1));
        ks{j}=alpha{j}.*0.20698;
        ksR(1,j) = mean(ks{j}(4500:end,1))./mean(ks{j}(100:500,1));
        lossT{j} = xL0.*gamma.^2.*beta{j}.*f0{j}./(alpha{j}.*f0{j});
        lossTR(1,j) = mean(lossT{j}(4500:end,1))./mean(lossT{j}(100:500,1));
        z{j} = 1:1:length(lossT{j});
    end
end
% function plotVar(Var1, Var2, varname1, varname2, plotname, fignum)
%         figure(fignum)
%         plot(Var1,Var2,'.', 'MarkerSize',5)
%         hold on
% %         %//number of points for the first part of the curve:
% %         n=300;
% %         %// Separate (x,y) into (x1,y1) and (x2,y2)
% %         Var11 = Var1(1:n); Var12=Var1(n+1:end);
% %         Var21 = Var2(1:n); Var22=Var2(n+1:end);
% % 
% %         %// fit a line y=A1*x+A2 to the first set of points:
% %         M=[Var11(:) ones(length(Var11),1)];
% %         A = M\Var21(:); %//A(1) is your slope, A(2) is your y-intercept
% % 
% %         %// fit a line y=B1*x+B2 to the second set of points:
% %         M=[Var12(:) ones(length(Var12),1)];
% %         B = M\Var22(:); %//B(1) is your slope, B(2) is your y-intercept
% % 
% %         %//Plot:
% %         plot(Var11,Var21,'o'); 
% %         plot (Var11, A(1)*Var11+A(2)) 
% %         plot(Var12,Var22,'o'); 
% %         plot (Var12, B(1)*Var12+B(2))
%         n=300;
%         F = @(B,Var1) min(B,Var1);
%         B_est = lsqcurvefit(F,0,Var1, Var2);
%         subplot(211)
%         plot(Var1,F(B_est,Var1),'r');
%         legend({'Original Data' char(F)});
%         title('Fitted function is not differentiable');
%         subplot(212)
%         ezplot(@(B) norm( F(B,Var1) - Var2 ).^2,[0 length(Var1)]), title('But error^2 is smooth');
%         hold on;
%         plot(B_est, norm( F(B_est,Var1) - Var2 ).^2,'ro')
% %         F = @(B,Var1) min(B(1),B(2)+B(3)*Var1);   %Form of the equation
% %             IC = [max(Var2) max(Var2) 0]; %Initial guess
% %             B = lsqcurvefit(F,IC,Var1,Var2,[min(Var2) -inf 0],[max(Var1) inf 0]);
% %             plot(Var1,F(B,Var1));
% %             a = (B(1) - B(2)) / B(3);
% %             cte = B(1);
% %             c = B(2);
% %             d = B(3);
%     hold off
%     title(plotname,'fontsize',15)
%     xlabel(varname1,'fontsize',15)
%     ylabel(varname2,'fontsize',15)
%     legend(gca,'1 mW','2 mW','4 mW','8 mW','12 mW')
%     set(gcf,'color','w');
%     %set(gca, 'YScale', 'log')
%     %set(gca, 'XScale', 'log')
%     %xlim([0 500])
% end

function plotVar(Var1, Var2, varname1, varname2, plotname,fignum)
    color={[8 150 255]./ 255,[240 120 0]./ 255,[251 200 0]./ 255, [156 0 145]./ 255,[10 150 10]./ 255};
    for j=1:length(Var1)
        trunc0=1;
        %Var1{j}=Var1{j}(trunc0:end);%52000
        %Var2{j}=Var2{j}(trunc0:end);
        %Var2{j}=Var2{j}-mean(Var2{j}(1:500))+mean(Var2{1}(1:500));
        figure(fignum)
        scatter(Var1{j},Var2{j},'.','MarkerEdgeAlpha',.1)
        hold on
    end    
    for i=1:length(Var1)
%         n = numel(Var1{i});
%         BIN_NUM = 30;
%         x_n = mean(reshape(Var1{i},BIN_NUM, n/BIN_NUM),1);
%         y_n = mean(reshape(Var2{i},BIN_NUM, n/BIN_NUM),1);
%         plot(x_n, y_n, 'k-','MarkerSize',1);
        M = movmean(Var2{i},[300 300]);
        trunc=500;
        alpha = .25;
        plot(Var1{i}(trunc:end),M(trunc:end),'Color',color{i},'LineStyle','-','LineWidth',1.8)%,'Color',[0,0,0]+alpha)
    end
    title(plotname,'fontsize',15)
    xlabel(varname1,'fontsize',15)
    ylabel(varname2,'fontsize',15)
    %legend(gca,'1 mW','2 mW','4 mW','8 mW','12 mW')
    %legend(gca,'0 \mum','1 \mum','5 \mum','10 \mum')
    set(gcf,'color','w');
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlim([10^(-2) 450])
    box on
    grid on
    hold off
end


