%Application 5: Battery state-of-charge and state-of-health estimation
close all
clc
clear

repeat=10000;%The number of times that the simulation is repeated at each measurement noise setup
% This simulation takes about 3 hours, you can reduce the repeat times to reduce the run time

%%
noise1_ref=0.001;
scale=-2:0.5:3;% The range of measurement noise is 10^scale
xplotEKF_1=[];
yplotEKF_1=[];
Time_EKF_1=[];
xplotEKF2_1=[];
yplotEKF2_1=[];
Time_EKF2_1=[];
xplotUKF_1=[];
yplotUKF_1=[];
Time_UKF_1=[];
xplotCKF_1=[];
yplotCKF_1=[];
Time_CKF_1=[];
xplotEKF_2=[];
yplotEKF_2=[];
Time_EKF_2=[];
xplotEKF2_2=[];
yplotEKF2_2=[];
Time_EKF2_2=[];
xplotUKF_2=[];
yplotUKF_2=[];
Time_UKF_2=[];
xplotCKF_2=[];
yplotCKF_2=[];
Time_CKF_2=[];
xplotIEKF_1=[];
yplotIEKF_1=[];
Time_IEKF_1=[];
xplotIEKF_2=[];
yplotIEKF_2=[];
Time_IEKF_2=[];
xplotIEKF2_1=[];
yplotIEKF2_1=[];
Time_IEKF2_1=[];
xplotIEKF2_2=[];
yplotIEKF2_2=[];
Time_IEKF2_2=[];
for i=scale
    magnitude=10^i*noise1_ref;
    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(0,magnitude,repeat,0);
    xplotEKF_1=[xplotEKF_1 SOC_rmse];
    yplotEKF_1=[yplotEKF_1 SOH_rmse];
    Time_EKF_1=[Time_EKF_1 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(1,magnitude,repeat,0);
    xplotEKF_2=[xplotEKF_2 SOC_rmse];
    yplotEKF_2=[yplotEKF_2 SOH_rmse];
    Time_EKF_2=[Time_EKF_2 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(0,magnitude,repeat,0.5);
    xplotIEKF_1=[xplotIEKF_1 SOC_rmse];
    yplotIEKF_1=[yplotIEKF_1 SOH_rmse];
    Time_IEKF_1=[Time_IEKF_1 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(1,magnitude,repeat,0.5);
    xplotIEKF_2=[xplotIEKF_2 SOC_rmse];
    yplotIEKF_2=[yplotIEKF_2 SOH_rmse];
    Time_IEKF_2=[Time_IEKF_2 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(0,magnitude,repeat,3);
    xplotEKF2_1=[xplotEKF2_1 SOC_rmse];
    yplotEKF2_1=[yplotEKF2_1 SOH_rmse];
    Time_EKF2_1=[Time_EKF2_1 Runtime];


    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(1,magnitude,repeat,3);
    xplotEKF2_2=[xplotEKF2_2 SOC_rmse];
    yplotEKF2_2=[yplotEKF2_2 SOH_rmse];
    Time_EKF2_2=[Time_EKF2_2 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(0,magnitude,repeat,3.5);
    xplotIEKF2_1=[xplotIEKF2_1 SOC_rmse];
    yplotIEKF2_1=[yplotIEKF2_1 SOH_rmse];
    Time_IEKF2_1=[Time_IEKF2_1 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(1,magnitude,repeat,3.5);
    xplotIEKF2_2=[xplotIEKF2_2 SOC_rmse];
    yplotIEKF2_2=[yplotIEKF2_2 SOH_rmse];
    Time_IEKF2_2=[Time_IEKF2_2 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(0,magnitude,repeat,1);
    xplotUKF_1=[xplotUKF_1 SOC_rmse];
    yplotUKF_1=[yplotUKF_1 SOH_rmse];
    Time_UKF_1=[Time_UKF_1 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(1,magnitude,repeat,1);
    xplotUKF_2=[xplotUKF_2 SOC_rmse];
    yplotUKF_2=[yplotUKF_2 SOH_rmse];
    Time_UKF_2=[Time_UKF_2 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(0,magnitude,repeat,2);
    xplotCKF_1=[xplotCKF_1 SOC_rmse];
    yplotCKF_1=[yplotCKF_1 SOH_rmse];
    Time_CKF_1=[Time_CKF_1 Runtime];

    [SOC_rmse,SOH_rmse,V_rmse, Runtime]=simulation(1,magnitude,repeat,2);
    xplotCKF_2=[xplotCKF_2 SOC_rmse];
    yplotCKF_2=[yplotCKF_2 SOH_rmse];
    Time_CKF_2=[Time_CKF_2 Runtime];
end
%% display Runtime
disp(['EKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_EKF_1)), ' ms'])
disp(['EKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_EKF_2)), ' ms'])
disp(['IEKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_IEKF_1)), ' ms'])
%disp(['IEKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_IEKF_2)), ' ms'])
disp(['EKF2 average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_EKF2_1)), ' ms'])
disp(['EKF2 average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_EKF2_2)), ' ms'])
%disp(['IEKF2 average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_IEKF2_1)), ' ms'])
%disp(['IEKF2 average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_IEKF2_2)), ' ms'])
disp(['UKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_UKF_1)), ' ms'])
disp(['UKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_UKF_2)), ' ms'])
disp(['CKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_CKF_1)), ' ms'])
disp(['CKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_CKF_2)), ' ms'])
%% figures

f1=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
h1_2=plot(noise1_ref*10.^scale, xplotEKF_2*100, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(noise1_ref*10.^scale, xplotEKF_1*100, '--o', 'DisplayName', "EKF (old)", 'Color', '#0072BD','LineWidth',1);
h5_1=plot(noise1_ref*10.^scale, xplotIEKF_1*100, '--d', 'DisplayName', "IEKF",'Color', '#4DBEEE','LineWidth',1);
h2_2=plot(noise1_ref*10.^scale, xplotEKF2_2*100, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(noise1_ref*10.^scale, xplotEKF2_1*100, '--s', 'DisplayName', "EKF2 (old)", 'Color', "#D95319",'LineWidth',1);
h3_2=plot(noise1_ref*10.^scale, xplotUKF_2*100, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(noise1_ref*10.^scale, xplotUKF_1*100, '--^', 'DisplayName', "UKF (old)", 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(noise1_ref*10.^scale, xplotCKF_2*100, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(noise1_ref*10.^scale, xplotCKF_1*100, '--x', 'DisplayName', "CKF (old)", 'Color', "#7E2F8E",'LineWidth',1);
%h6_2=plot(noise1_ref*10.^scale, xplotIEKF2_2*100, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(noise1_ref*10.^scale, xplotIEKF2_1*100, '--v', 'Color', "#A2142F",'LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
%xlabel('Measurement noise (per unit)',FontSize=13)
ylabel('SOC RMSE (%)',FontSize=13)
grid on
ylim([0.001 10])
xlim([10^-5 1])
set(gca, 'YTick', [0.001 0.01 0.1 1 10]);
%legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast',FontSize=11)

nexttile
hold on
h1_2=plot(noise1_ref*10.^scale, yplotEKF_2*100, '-o', 'DisplayName', "EKF (new)", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(noise1_ref*10.^scale, yplotEKF_1*100, '--o', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(noise1_ref*10.^scale, yplotEKF2_2*100, '-s', 'DisplayName', "EKF2 (new)", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(noise1_ref*10.^scale, yplotEKF2_1*100, '--s', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(noise1_ref*10.^scale, yplotUKF_2*100, '-^', 'DisplayName', "UKF (new)", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(noise1_ref*10.^scale, yplotUKF_1*100, '--^', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(noise1_ref*10.^scale, yplotCKF_2*100, '-x', 'DisplayName', "CKF (new)", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(noise1_ref*10.^scale, yplotCKF_1*100, '--x', 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(noise1_ref*10.^scale, yplotIEKF_2*100, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
h5_1=plot(noise1_ref*10.^scale, yplotIEKF_1*100, '--d', 'Color', '#4DBEEE','LineWidth',1);
%h6_2=plot(noise1_ref*10.^scale, yplotIEKF2_2*100, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(noise1_ref*10.^scale, yplotIEKF2_1*100, '--v', 'Color', "#A2142F",'LineWidth',1);
h1_1=plot(nan, nan, 'o', 'Color', '#0072BD', 'DisplayName', 'EKF','LineWidth',1);
h2_1=plot(nan, nan, 's', 'Color', "#D95319", 'DisplayName', 'EKF2','LineWidth',1);
h3_1=plot(nan, nan, '^', 'Color', "#EDB120", 'DisplayName', 'UKF','LineWidth',1);
h4_1=plot(nan, nan, 'x', 'Color', "#7E2F8E", 'DisplayName', 'CKF','LineWidth',1);
h5_1=plot(nan, nan, 'd', 'Color', "#4DBEEE", 'DisplayName', 'IEKF','LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
xlabel('Measurement standard deviation (V)',FontSize=13)
ylabel('SOH RMSE (%)',FontSize=13)
legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast',FontSize=11)
grid on
ylim([0.1 12])
set(gca, 'YTick', [0.01 0.1 1 10]);
xlim([10^-5 1])
% exportgraphics(f1,'SOC-SOH.png','Resolution',900)
%%
% figure
% f2=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
% nexttile
% hold on
% h1_2=plot(noise1_ref*10.^scale, xplotEKF_2*100, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(noise1_ref*10.^scale, xplotEKF_1*100, '--o', 'Color', '#0072BD','LineWidth',1);
% h5_2=plot(noise1_ref*10.^scale, xplotIEKF_2*100, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
% h5_1=plot(noise1_ref*10.^scale, xplotIEKF_1*100, '--d', 'Color', '#4DBEEE','LineWidth',1);
% h2_2=plot(noise1_ref*10.^scale, xplotEKF2_2*100, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(noise1_ref*10.^scale, xplotEKF2_1*100, '--s', 'Color', "#D95319",'LineWidth',1);
% h6_2=plot(noise1_ref*10.^scale, xplotIEKF2_2*100, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(noise1_ref*10.^scale, xplotIEKF2_1*100, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
% %xlabel('Measurement noise (per unit)',FontSize=13)
% ylabel('SOC RMSE (%)',FontSize=13)
% grid on
% ylim([0.001 10])
% xlim([10^-5 1])
% set(gca, 'YTick', [0.001 0.01 0.1 1 10]);
% %legend([h1_2 h2_2 h3_2 h4_2], 'Location','southeast',FontSize=10)
% 
% nexttile
% hold on
% h1_2=plot(noise1_ref*10.^scale, yplotEKF_2*100, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(noise1_ref*10.^scale, yplotEKF_1*100, '--o', 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(noise1_ref*10.^scale, yplotEKF2_2*100, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(noise1_ref*10.^scale, yplotEKF2_1*100, '--s', 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(noise1_ref*10.^scale, yplotIEKF2_2*100, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
% h5_1=plot(noise1_ref*10.^scale, yplotIEKF2_1*100, '--d', 'Color', '#4DBEEE','LineWidth',1);
% h6_2=plot(noise1_ref*10.^scale, yplotIEKF2_2*100, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(noise1_ref*10.^scale, yplotIEKF2_1*100, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
% xlabel('Measurement standard deviation (V)',FontSize=13)
% ylabel('SOH RMSE (%)',FontSize=13)
% legend([h1_2 h2_2 h5_2 h6_2], 'Location','southeast',FontSize=12)
% grid on
% ylim([0.1 12])
% set(gca, 'YTick', [0.01 0.1 1 10]);
% xlim([10^-5 1])
% exportgraphics(f2,'SOC-SOH-IKF.png','Resolution',900)
function [SOC_rmse,SOH_rmse,V_rmse,Runtime]=simulation(improvement1,noise,repeat,KFtype)
%% Problem setup
%a_100 are the parameters at 100% SOH
a_100=[1390.38022324692	-6961.31052245824	14760.3117954132	-17230.9243091027	12055.7081234694	-5162.75566125427	1330.60562294128	-196.376699215906	15.5992721935811	2.96139977813363];
%a_80 are the parameters at 80% SOH
a_80=[813.942183091546	-4229.96630149708	9345.48940409742	-11415.3796212382	8396.14648778192	-3801.06696008866	1043.08567316897	-165.294460011886	14.2780123166013	2.96071001192947];
%OCV = a0 + a1*SOC + a2*SOC^2 + a3*SOC^3 + ... + a9*SOC^9
% a0,...,a9 are all linear functions of SOH

delta_t=1;
T1=15;%initial rest
T2=60;%discharge
T3=30;%rest
T4=60;%charge
T5=15;%rest
I0=2;% 2C charge & discharge
R1=0.08;
R2=0.05;
R2C_reciporal=0.008;
I_error=0.001;
Tlimit=T1+T2+T3+T4+T5;
SOC_rmse=0;
SOH_rmse=0;
V_rmse=0;
rng(123)
Runtime=0;
for iii=1:repeat
    Q0=3600; %1Ah
    states_true=[0.6;0;0.9];% SOC V SOH
    states_est=[0.8;0;1];
    Variance=[0.2^2 0 0; 0 1e-5^2 0; 0 0 0.1^2];
    %% generate profile
    %current profile: 15s rest, 60s discharge, 30s rest, 60s charge, 15s rest
    Measurementnoise=noise^2;
    I_true=0;
    t=delta_t:delta_t:Tlimit;
    stepnum=length(t);
    for i=1:stepnum
        states=states_true(:,i);
        if(delta_t*i<=T1)
            I=0;
        elseif(delta_t*i<=T1+T2)
            I=-I0;
        elseif(delta_t*i<=T1+T2+T3)
            I=0;
        elseif(delta_t*i<=T1+T2+T3+T4)
            I=I0;
        else
            I=0;
        end
        states=[states(1)+I*delta_t/states(3)/Q0;
            states(2)*exp(-delta_t*R2C_reciporal)+(R2-R2*exp(-delta_t*R2C_reciporal))*I;
            states(3)];
        states_true=[states_true states];
        I_true=[I_true I];
    end
    a_true=a_100*(states_true(3)-0.8)*5+a_80*(1-states_true(3))*5;
    Measurements=OCV(states_true(1,:),a_true)+states_true(2,:)+I_true*R1;
    Measurements_noisy=Measurements+randn(1,stepnum+1)*noise;
    I_noisy=I_true+randn(1,stepnum+1)*I_error;
    %% KF
    tic;
    for i=1:stepnum
        if(KFtype==0)
            states=states_est(:,i);
            F=[1 0 -delta_t*I_noisy(i+1)/states(3)^2/Q0; 0 exp(-delta_t*R2C_reciporal) 0; 0 0 1];
            states=[states(1)+I_noisy(i+1)*delta_t/states(3)/Q0;
            states(2)*exp(-delta_t*R2C_reciporal)+(R2-R2*exp(-delta_t*R2C_reciporal))*I_noisy(i+1);
            states(3)];
            B = [delta_t/states(3)/Q0; R2-R2*exp(-delta_t*R2C_reciporal); 0];
            Q=B*B.'*I_error^2;
            Variance=F*Variance*F.'+Q;
            a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
            residual=Measurements_noisy(:,i+1)-OCV(states(1),a_est)-states(2)-I_noisy(i+1)*R1;
            H=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
            S=H*Variance*H.'+Measurementnoise+R1^2*I_error^2;
            K=Variance*H.'*S^(-1);
            states0=states;
            states=states+K*residual;
            if(improvement1==1)
                a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
                H2=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
                temp=Variance-K*H2*Variance-Variance*H2'*K'+K*(H2*Variance*H2'+Measurementnoise+R1^2*I_error^2)*K';
                if(sum(diag(temp)))<sum(diag(Variance))
                    Variance=temp;
                else
                    states=states0;
                end
            else
                Variance=Variance-K*H*Variance;
            end
            states_est(:,i+1)=states;
        elseif(KFtype==0.5)
            states=states_est(:,i);
            F=[1 0 -delta_t*I_noisy(i+1)/states(3)^2/Q0; 0 exp(-delta_t*R2C_reciporal) 0; 0 0 1];
            states=[states(1)+I_noisy(i+1)*delta_t/states(3)/Q0;
            states(2)*exp(-delta_t*R2C_reciporal)+(R2-R2*exp(-delta_t*R2C_reciporal))*I_noisy(i+1);
            states(3)];
            B = [delta_t/states(3)/Q0; R2-R2*exp(-delta_t*R2C_reciporal); 0];
            Q=B*B.'*I_error^2;
            Variance=F*Variance*F.'+Q;
            a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
            residual=Measurements_noisy(:,i+1)-OCV(states(1),a_est)-states(2)-I_noisy(i+1)*R1;
            H=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
            S=H*Variance*H.'+Measurementnoise+R1^2*I_error^2;
            K=Variance*H.'*S^(-1);
            states0=states;
            states=states+K*residual;
            steplennow = norm(K*residual);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
                residual = Measurements_noisy(:,i+1)-OCV(states(1),a_est)-states(2)-I_noisy(i+1)*R1;
                H2=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
                S2 = H2*Variance*H2.'+Measurementnoise+R1^2*I_error^2;
                K2 = Variance*H2.'*S2^(-1);
                dx = states0 + K2*(residual-H2*(states0-states))-states;
                steplen_previous=steplennow;
                steplennow=norm(dx);
                if(steplen_previous<steplennow)
                    break;
                else
                    change = max(abs(dx./states));
                    states = states + dx;
                    K = K2;
                    H = H2;
                end
            end
            if(improvement1==1)
                a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
                H2=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
                temp=Variance-K*H2*Variance-Variance*H2'*K'+K*(H2*Variance*H2'+Measurementnoise+R1^2*I_error^2)*K';
                if(sum(diag(temp))<sum(diag(Variance)))
                    Variance=temp;
                else
                    states=states0;
                end
            else
                Variance=Variance-K*H*Variance;
            end
            states_est(:,i+1)=states;
        elseif(KFtype==1)
            L=real(sqrtm(Variance));
            state=states_est(:,i);
            n=3;
            lambda=(1e-6-1)*n;
            alpha=1e-3;
            beta=2;
            mnum=1;
            states=zeros(n,n*2+1);
            states(:,1)=state;
            for ii=1:n
                states(:,1+ii)=state+sqrt(lambda+n)*L(:,ii);
                states(:,n+ii+1)=state-sqrt(lambda+n)*L(:,ii);  
            end
            for ii=1:2*n+1
                states(:,ii)=[states(1,ii)+I_noisy(i+1)*delta_t/states(3,ii)/Q0;
            states(2,ii)*exp(-delta_t*R2C_reciporal)+(R2-R2*exp(-delta_t*R2C_reciporal))*I_noisy(i+1);
            states(3,ii)];
                if(ii==1)
                    state=states(:,ii)*lambda/(n+lambda);
                else
                    state=state+states(:,ii)/(n+lambda)/2;
                end
            end
            B = [delta_t/state(3)/Q0; R2-R2*exp(-delta_t*R2C_reciporal); 0];
            Q=B*B.'*I_error^2;
            Variance=(states(:,1)-state)*(states(:,1)-state).'*(lambda/(n+lambda)+1-alpha^2+beta)+Q;
            for ii=2:2*n+1
                Variance=Variance+(state-states(:,ii))*(state-states(:,ii)).'/(n+lambda)/2;
            end
            % obtain measurement
            z = Measurements_noisy(:,i+1);
            % Predict Measurement From Propagated Sigma Points
            measures=zeros(mnum,2*n+1);
            for ii=1:2*n+1
                a_est=a_100*(states(3,ii)-0.8)*5+a_80*(1-states(3,ii))*5;
                measures(:,ii) = OCV(states(1,ii),a_est)+states(2,ii)+I_noisy(i+1)*R1; % predicted measurement
                if(ii==1)
                    m_exp=lambda/(n+lambda)*measures(:,ii);
                else
                    m_exp=m_exp+1/(n+lambda)/2*measures(:,ii);
                end
            end
            % Estimate Mean And Covariance of Predicted Measurement
            Py=(lambda/(n+lambda)+1-alpha^2+beta)*(measures(:,1)-m_exp)*(measures(:,1)-m_exp).'+Measurementnoise+R1^2*I_error^2;
            Pxy=(lambda/(n+lambda)+1-alpha^2+beta)*(states(:,1)-state)*(measures(:,1)-m_exp).';
            for ii=2:2*n+1
                Py=Py+1/(n+lambda)/2*(measures(:,ii)-m_exp)*(measures(:,ii)-m_exp).';
                Pxy=Pxy+1/(n+lambda)/2*(states(:,ii)-state)*(measures(:,ii)-m_exp).';
            end
            %kalman gain
            K=Pxy/Py;
            %update
            state0=state;
            state=state+K*(z-m_exp);
            if(improvement1==1)
                L=real(sqrtm(Variance));
                states=zeros(n,n*2+1);
                states(:,1)=state;
                for ii=1:n
                    states(:,1+ii)=state+sqrt(lambda+n)*L(:,ii);
                    states(:,n+ii+1)=state-sqrt(lambda+n)*L(:,ii);
                end
                measures=zeros(mnum,2*n+1);
                for ii=1:2*n+1
                    a_est=a_100*(states(3,ii)-0.8)*5+a_80*(1-states(3,ii))*5;
                    measures(:,ii) = OCV(states(1,ii),a_est)+states(2,ii)+I_noisy(i+1)*R1; % predicted measurement
                    if(ii==1)
                        m_exp=lambda/(n+lambda)*measures(:,ii);
                    else
                        m_exp=m_exp+1/(n+lambda)/2*measures(:,ii);
                    end
                end
                % Estimate Mean And Covariance of Predicted Measurement
                Py=(lambda/(n+lambda)+1-alpha^2+beta)*(measures(:,1)-m_exp)*(measures(:,1)-m_exp).'+Measurementnoise+R1^2*I_error^2;
                Pxy=(lambda/(n+lambda)+1-alpha^2+beta)*(states(:,1)-state)*(measures(:,1)-m_exp).';
                for ii=2:2*n+1
                    Py=Py+1/(n+lambda)/2*(measures(:,ii)-m_exp)*(measures(:,ii)-m_exp).';
                    Pxy=Pxy+1/(n+lambda)/2*(states(:,ii)-state)*(measures(:,ii)-m_exp).';
                end
                temp=Variance+K*Py*K.'-Pxy*K.'-K*Pxy.';
                if(sum(diag(temp))<sum(diag(Variance)))
                    Variance=temp;
                else
                    state=state0;
                end
            else
                Variance=Variance-K*Py*K.';
            end
            states_est(:,i+1)=state;   
        elseif(KFtype==2)
            %CKF
            L=real(sqrtm(Variance));
            state=states_est(:,i);
            n=3;
            mnum=1;
            states=zeros(n,n*2);
            for ii=1:n
                states(:,ii)=state+sqrt(n)*L(:,ii);
                states(:,n+ii)=state-sqrt(n)*L(:,ii);  
            end
            state=0;
            for ii=1:2*n
                states(:,ii)=[states(1,ii)+I_noisy(i+1)*delta_t/states(3,ii)/Q0;
            states(2,ii)*exp(-delta_t*R2C_reciporal)+(R2-R2*exp(-delta_t*R2C_reciporal))*I_noisy(i+1);
            states(3,ii)];
                state=state+states(:,ii)/n/2;
            end
            B = [delta_t/state(3)/Q0; R2-R2*exp(-delta_t*R2C_reciporal); 0];
            Q=B*B.'*I_error^2;
            Variance=Q;
            for ii=1:2*n
                Variance=Variance+(state-states(:,ii))*(state-states(:,ii)).'/n/2;
            end
            % obtain measurement
            z = Measurements_noisy(:,i+1);
            % Predict Measurement From Propagated Sigma Points
            L=real(sqrtm(Variance));
            states=zeros(n,n*2);
            for ii=1:n
                states(:,ii)=state+sqrt(n)*L(:,ii);
                states(:,n+ii)=state-sqrt(n)*L(:,ii);
            end
            measures=zeros(mnum,2*n);
            m_exp=0;
            for ii=1:2*n
                a_est=a_100*(states(3,ii)-0.8)*5+a_80*(1-states(3,ii))*5;
                measures(:,ii) = OCV(states(1,ii),a_est)+states(2,ii)+I_noisy(i+1)*R1; % predicted measurement
                m_exp=m_exp+1/n/2*measures(:,ii);
            end
            % Estimate Mean And Covariance of Predicted Measurement
            Py=Measurementnoise+R1^2*I_error^2;
            Pxy=0;
            for ii=1:2*n
                Py=Py+1/n/2*(measures(:,ii)-m_exp)*(measures(:,ii)-m_exp).';
                Pxy=Pxy+1/n/2*(states(:,ii)-state)*(measures(:,ii)-m_exp).';
            end
            %kalman gain
            K=Pxy/Py;
            %update
            state0=state;
            state=state+K*(z-m_exp);
            if(improvement1==1)
                states=zeros(n,n*2);
                for ii=1:n
                    states(:,ii)=state+sqrt(n)*L(:,ii);
                    states(:,n+ii)=state-sqrt(n)*L(:,ii);
                end
                measures=zeros(mnum,2*n);
                m_exp=0;
                for ii=1:2*n
                    a_est=a_100*(states(3,ii)-0.8)*5+a_80*(1-states(3,ii))*5;
                    measures(:,ii) = OCV(states(1,ii),a_est)+states(2,ii)+I_noisy(i+1)*R1; % predicted measurement
                    m_exp=m_exp+1/n/2*measures(:,ii);
                end
                % Estimate Mean And Covariance of Predicted Measurement
                Py=Measurementnoise+R1^2*I_error^2;
                Pxy=0;
                for ii=1:2*n
                    Py=Py+1/n/2*(measures(:,ii)-m_exp)*(measures(:,ii)-m_exp).';
                    Pxy=Pxy+1/n/2*(states(:,ii)-state)*(measures(:,ii)-m_exp).';
                end
                temp=Variance+K*Py*K.'-Pxy*K.'-K*Pxy.';
                if(sum(diag(temp))<sum(diag(Variance)))
                    Variance=temp;
                else
                    state=state0;
                end
            else
                Variance=Variance-K*Py*K.';
            end
            states_est(:,i+1)=state;   
        elseif(KFtype==3)
            %2-EKF
            states=states_est(:,i);
            states=[states(1)+I_noisy(i+1)*delta_t/states(3)/Q0;
            states(2)*exp(-delta_t*R2C_reciporal)+(R2-R2*exp(-delta_t*R2C_reciporal))*I_noisy(i+1);
            states(3)];
            n=length(states);
            Fxx=zeros(n,n,n);
            deltaP=zeros(n,n);
            Fxx(3,3,1)=2*I_noisy(i+1)*delta_t/states(3)^3/Q0;
            for ii=1:n
                states(ii)=states(ii)+0.5*sum(diag(Fxx(:,:,ii)*Variance));
                for jj=1:n
                    deltaP(ii,jj)=sum(diag(Fxx(:,:,ii)*Variance*Fxx(:,:,jj)*Variance));
                end
            end
            F=[1 0 -delta_t*I_noisy(i+1)/states(3)^2/Q0; 0 exp(-delta_t*R2C_reciporal) 0; 0 0 1];
            B = [delta_t/states(3)/Q0; R2-R2*exp(-delta_t*R2C_reciporal); 0];
            Q=B*B.'*I_error^2;
            Variance=F*Variance*F.'+0.5*deltaP+Q;
            %update
            mnum=1;
            a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=[d2OCVd2SOC(states(1),a_est) 0 dOCVdSOC(states(1),a_100-a_80)*5;
                0 0 0;
                dOCVdSOC(states(1),a_100-a_80)*5 0 0];
            H=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
            S=H*Variance*H.'+Measurementnoise+R1^2*I_error^2;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*Variance*Hxx(:,:,jj)*Variance));
                end
            end
            K=Variance*H.'*S^(-1);
            residual=Measurements_noisy(:,i+1)-OCV(states(1),a_est)-states(2)-I_noisy(i+1)*R1;
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));
            end
            states0=states;
            states=states+K*residual;
            if(improvement1==1)
                a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[d2OCVd2SOC(states(1),a_est) 0 dOCVdSOC(states(1),a_100-a_80)*5;
                    0 0 0;
                    dOCVdSOC(states(1),a_100-a_80)*5 0 0];
                H2=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
                S2=H2*Variance*H2.'+Measurementnoise+R1^2*I_error^2;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*Variance*Hxx2(:,:,jj)*Variance));
                    end
                end
                temp=Variance+K*S2*K.'-Variance*H2.'*K.'-K*H2*Variance;
                if(sum(diag(temp))<sum(diag(Variance)))
                    Variance=temp;
                else
                    states=states0;
                end
            else
                Variance=Variance-K*S*K.';
            end
            states_est(:,i+1)=states;
        else
            %2-IEKF
            states=states_est(:,i);
            states=[states(1)+I_noisy(i+1)*delta_t/states(3)/Q0;
            states(2)*exp(-delta_t*R2C_reciporal)+(R2-R2*exp(-delta_t*R2C_reciporal))*I_noisy(i+1);
            states(3)];
            n=length(states);
            Fxx=zeros(n,n,n);
            deltaP=zeros(n,n);
            Fxx(3,3,1)=2*I_noisy(i+1)*delta_t/states(3)^3/Q0;
            for ii=1:n
                states(ii)=states(ii)+0.5*sum(diag(Fxx(:,:,ii)*Variance));
                for jj=1:n
                    deltaP(ii,jj)=sum(diag(Fxx(:,:,ii)*Variance*Fxx(:,:,jj)*Variance));
                end
            end
            F=[1 0 -delta_t*I_noisy(i+1)/states(3)^2/Q0; 0 exp(-delta_t*R2C_reciporal) 0; 0 0 1];
            B = [delta_t/states(3)/Q0; R2-R2*exp(-delta_t*R2C_reciporal); 0];
            Q=B*B.'*I_error^2;
            Variance=F*Variance*F.'+0.5*deltaP+Q;
            %update
            mnum=1;
            a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=[d2OCVd2SOC(states(1),a_est) 0 dOCVdSOC(states(1),a_100-a_80)*5;
                0 0 0;
                dOCVdSOC(states(1),a_100-a_80)*5 0 0];
            H=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
            S=H*Variance*H.'+Measurementnoise+R1^2*I_error^2;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*Variance*Hxx(:,:,jj)*Variance));
                end
            end
            K=Variance*H.'*S^(-1);
            residual=Measurements_noisy(:,i+1)-OCV(states(1),a_est)-states(2)-I_noisy(i+1)*R1;
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));
            end
            states0=states;
            states=states+K*residual;
            
            steplennow = norm(K*residual);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
                residual = Measurements_noisy(:,i+1)-OCV(states(1),a_est)-states(2)-I_noisy(i+1)*R1;
                H2=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[d2OCVd2SOC(states(1),a_est) 0 dOCVdSOC(states(1),a_100-a_80)*5;
                0 0 0;
                dOCVdSOC(states(1),a_100-a_80)*5 0 0];
                for ii=1:mnum
                    residual(ii)=residual(ii)-0.5*(states0-states).'*Hxx2(:,:,ii)*(states0-states);
                    residual(ii)=residual(ii)-0.5*sum(diag(Hxx2(:,:,ii)*Variance));
                end
                S2 = H2*Variance*H2.'+Measurementnoise+R1^2*I_error^2;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*Variance*Hxx2(:,:,jj)*Variance));
                    end
                end
                K2 = Variance*H2.'*S2^(-1);
                dx = states0 + K2*(residual-H2*(states0-states))-states;
                steplen_previous=steplennow;
                steplennow=norm(dx);
                if(steplen_previous<steplennow)
                    break;
                else
                    change = max(abs(dx./states));
                    states = states + dx;
                    K = K2;
                    S = S2;
                end
            end
            if(improvement1==1)
                a_est=a_100*(states(3)-0.8)*5+a_80*(1-states(3))*5;
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[d2OCVd2SOC(states(1),a_est) 0 dOCVdSOC(states(1),a_100-a_80)*5;
                    0 0 0;
                    dOCVdSOC(states(1),a_100-a_80)*5 0 0];
                H2=[dOCVdSOC(states(1),a_est) 1 OCV(states(1),a_100-a_80)*5];
                S2=H2*Variance*H2.'+Measurementnoise+R1^2*I_error^2;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*Variance*Hxx2(:,:,jj)*Variance));
                    end
                end
                temp=Variance+K*S2*K.'-Variance*H2.'*K.'-K*H2*Variance;
                if(sum(diag(temp))<sum(diag(Variance)))
                    Variance=temp;
                else
                    states=states0;
                end
            else
                Variance=Variance-K*S*K.';
            end
            states_est(:,i+1)=states;
        end
    end
    temp = toc;
    Runtime = Runtime + temp/repeat;
    SOC_est=states_est(1,:);
    V_est=states_est(2,:);
    SOH_est=states_est(3,:);
    %x_rmse=x_rmse+rmse(x,x_est)/repeat;
    %y_rmse=y_rmse+rmse(y,y_est)/repeat;
    %phi_rmse=phi_rmse+180/pi*rmse(states_true(3,:),states_est(3,:))/repeat;
    SOC_rmse=SOC_rmse+(states_true(1,length(SOC_est))-SOC_est(length(SOC_est)))^2/repeat;
    V_rmse=V_rmse+(states_true(2,length(SOC_est))-V_est(length(SOC_est)))^2/repeat;
    SOH_rmse=SOH_rmse+(states_true(3,length(SOC_est))-SOH_est(length(SOC_est)))^2/repeat;   
end
SOC_rmse=sqrt(SOC_rmse);
V_rmse=sqrt(V_rmse);
SOH_rmse=sqrt(SOH_rmse);
end

function dVdSOC = dOCVdSOC(SOC,coeficients)
xspace=zeros(length(coeficients),length(SOC));
for i=1:length(coeficients)-1
    xspace(i,:)=(length(coeficients)-i).*SOC.^(length(coeficients)-i-1);
end
dVdSOC=coeficients*xspace;
end

function d2Vd2SOC = d2OCVd2SOC(SOC,coeficients)
xspace=zeros(length(coeficients),length(SOC));
for i=1:length(coeficients)-2
    xspace(i,:)=(length(coeficients)-i).*(length(coeficients)-i-1).*SOC.^(length(coeficients)-i-2);
end
d2Vd2SOC=coeficients*xspace;
end

function V = OCV(SOC,coeficients)
xspace=NaN(length(coeficients),length(SOC));
for i=1:length(coeficients)
    xspace(i,:)=SOC.^(length(coeficients)-i);
end
V=coeficients*xspace;
end
