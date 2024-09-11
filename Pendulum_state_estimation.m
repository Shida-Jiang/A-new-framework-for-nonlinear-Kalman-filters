%Application 4: pendulum state estimation

close all
clc
clear
repeat=10000; %The number of times that the simulation is repeated at each measurement noise setup
% This simulation takes about 30 minutes, you can reduce the repeat times to reduce the run time
%%
noise1_ref=1;
scale=-4:0.5:1; % The range of measurement noise is 10^scale
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
    noise1=noise1_ref*10^i;

    [x_rmse,y_rmse,Runtime]=simulation(0,noise1,repeat,0.5);
    xplotIEKF_1=[xplotIEKF_1 x_rmse];
    yplotIEKF_1=[yplotIEKF_1 y_rmse];
    Time_IEKF_1=[Time_IEKF_1 Runtime];
    [x_rmse,y_rmse,Runtime]=simulation(1,noise1,repeat,0.5);
    xplotIEKF_2=[xplotIEKF_2 x_rmse];
    yplotIEKF_2=[yplotIEKF_2 y_rmse];
    Time_IEKF_2=[Time_IEKF_2 Runtime];

    [x_rmse,y_rmse,Runtime]=simulation(0,noise1,repeat,0);
    xplotEKF_1=[xplotEKF_1 x_rmse];
    yplotEKF_1=[yplotEKF_1 y_rmse];
    Time_EKF_1=[Time_EKF_1 Runtime];
    [x_rmse,y_rmse,Runtime]=simulation(1,noise1,repeat,0);
    xplotEKF_2=[xplotEKF_2 x_rmse];
    yplotEKF_2=[yplotEKF_2 y_rmse];
    Time_EKF_2=[Time_EKF_2 Runtime];

    [x_rmse,y_rmse,Runtime]=simulation(0,noise1,repeat,3);
    xplotEKF2_1=[xplotEKF2_1 x_rmse];
    yplotEKF2_1=[yplotEKF2_1 y_rmse];
    Time_EKF2_1=[Time_EKF2_1 Runtime];
    [x_rmse,y_rmse,Runtime]=simulation(1,noise1,repeat,3);
    xplotEKF2_2=[xplotEKF2_2 x_rmse];
    yplotEKF2_2=[yplotEKF2_2 y_rmse];
    Time_EKF2_2=[Time_EKF2_2 Runtime];

    [x_rmse,y_rmse,Runtime]=simulation(0,noise1,repeat,3.5);
    xplotIEKF2_1=[xplotIEKF2_1 x_rmse];
    yplotIEKF2_1=[yplotIEKF2_1 y_rmse];
    Time_IEKF2_1=[Time_IEKF2_1 Runtime];
    [x_rmse,y_rmse,Runtime]=simulation(1,noise1,repeat,3.5);
    xplotIEKF2_2=[xplotIEKF2_2 x_rmse];
    yplotIEKF2_2=[yplotIEKF2_2 y_rmse];
    Time_IEKF2_2=[Time_IEKF2_2 Runtime];

    [x_rmse,y_rmse,Runtime]=simulation(0,noise1,repeat,1);
    xplotUKF_1=[xplotUKF_1 x_rmse];
    yplotUKF_1=[yplotUKF_1 y_rmse];
    Time_UKF_1=[Time_UKF_1 Runtime];
    [x_rmse,y_rmse,Runtime]=simulation(1,noise1,repeat,1);
    xplotUKF_2=[xplotUKF_2 x_rmse];
    yplotUKF_2=[yplotUKF_2 y_rmse];
    Time_UKF_2=[Time_UKF_2 Runtime];

    [x_rmse,y_rmse,Runtime]=simulation(0,noise1,repeat,2);
    xplotCKF_1=[xplotCKF_1 x_rmse];
    yplotCKF_1=[yplotCKF_1 y_rmse];
    Time_CKF_1=[Time_CKF_1 Runtime];
    [x_rmse,y_rmse,Runtime]=simulation(1,noise1,repeat,2);
    xplotCKF_2=[xplotCKF_2 x_rmse];
    yplotCKF_2=[yplotCKF_2 y_rmse];
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
h1_2=plot(10.^scale, yplotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, yplotEKF_1, '--o', 'DisplayName', "EKF (old)", 'Color', '#0072BD','LineWidth',1);
%h5_2=plot(10.^scale, yplotIEKF_2, '-d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
h5_1=plot(10.^scale, yplotIEKF_1, '--d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
h2_2=plot(10.^scale, yplotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, yplotEKF2_1, '--s', 'DisplayName', "EKF2 (old)", 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, yplotUKF_2, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, yplotUKF_1, '--^', 'DisplayName', "UKF (old)", 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, yplotCKF_2, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, yplotCKF_1, '--x', 'DisplayName', "CKF (old)", 'Color', "#7E2F8E",'LineWidth',1);
%h6_2=plot(10.^scale, yplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, yplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
%xlabel('Measurement noise (per unit)',FontSize=12)
ylabel('Angle RMSE (rad)',FontSize=12)
h1_1=plot(nan, nan, '--', 'Color', 'k', 'DisplayName', 'Old framework','LineWidth',1);
h2_1=plot(nan, nan, 'Color', "k", 'DisplayName', 'New framework','LineWidth',1);
leg=legend([h1_1 h2_1], 'Location','southeast',FontSize=11);
title(leg,'Line styles')
xlim([10^-4 10^1])
set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
ylim([1e-6 1])
set(gca, 'YTick', [1e-6 1e-4 1e-2 1]);
grid on
nexttile
hold on
h1_2=plot(10.^scale, xplotEKF_2, '-o', 'DisplayName', "EKF (new)", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, xplotEKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, xplotEKF2_2, '-s', 'DisplayName', "EKF2 (new)", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, xplotEKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, xplotUKF_2, '-^', 'DisplayName', "UKF (new)", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, xplotUKF_1, '--^', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, xplotCKF_2, '-x', 'DisplayName', "CKF (new)", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, xplotCKF_1, '--x', 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, xplotIEKF_2, '-d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
h5_1=plot(10.^scale, xplotIEKF_1, '--d', 'Color', "#4DBEEE",'LineWidth',1);
%h6_2=plot(10.^scale, xplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, xplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
h1_1=plot(nan, nan, 'o', 'Color', '#0072BD', 'DisplayName', 'EKF','LineWidth',1);
h2_1=plot(nan, nan, 's', 'Color', "#D95319", 'DisplayName', 'EKF2','LineWidth',1);
h3_1=plot(nan, nan, '^', 'Color', "#EDB120", 'DisplayName', 'UKF','LineWidth',1);
h4_1=plot(nan, nan, 'x', 'Color', "#7E2F8E", 'DisplayName', 'CKF','LineWidth',1);
h5_1=plot(nan, nan, 'd', 'Color', "#4DBEEE", 'DisplayName', 'IEKF','LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
xlabel('Measurement standard deviation (N)',FontSize=12)
ylabel('Anglular Speed RMSE (rad/s)',FontSize=12)
leg=legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast','NumColumns',2,FontSize=11);
title(leg,'Line colors')
xlim([10^-4 10^1])
set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
ylim([1e-6 1])
set(gca, 'YTick', [1e-6 1e-4 1e-2 1]);
grid on
%exportgraphics(f1,'pendulum.png','Resolution',900)
%%
% figure
% f2=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
% nexttile
% hold on
% h1_2=plot(10.^scale, yplotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, yplotEKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
% h5_2=plot(10.^scale, yplotIEKF_2, '-d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
% h5_1=plot(10.^scale, yplotIEKF_1, '--d', 'Color', "#4DBEEE",'LineWidth',1);
% h2_2=plot(10.^scale, yplotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, yplotEKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
% h6_2=plot(10.^scale, yplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, yplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
% %xlabel('Measurement noise (per unit)',FontSize=12)
% ylabel('Angle RMSE (rad)',FontSize=12)
% %legend([h1_2 h2_2 h3_2 h4_2], 'Location','southeast',FontSize=10)
% xlim([10^-4 10^1])
% set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
% ylim([1e-6 1])
% set(gca, 'YTick', [1e-6 1e-4 1e-2 1]);
% grid on
% nexttile
% hold on
% h1_2=plot(10.^scale, xplotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, xplotEKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(10.^scale, xplotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, xplotEKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(10.^scale, xplotIEKF_2, '-d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
% h5_1=plot(10.^scale, xplotIEKF_1, '--d', 'Color', "#4DBEEE",'LineWidth',1);
% h6_2=plot(10.^scale, xplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, xplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
% xlabel('Measurement standard deviation (N)',FontSize=12)
% ylabel('Anglular Speed RMSE (rad/s)',FontSize=12)
% legend([h1_2 h2_2 h5_2 h6_2], 'Location','southeast',FontSize=12)
% xlim([10^-4 10^1])
% set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
% ylim([1e-6 1])
% set(gca, 'YTick', [1e-6 1e-4 1e-2 1]);
% grid on
% exportgraphics(f2,'pendulum-IKF.png','Resolution',900)
function [sitadot_rmse,sita_rmse,Runtime]=simulation(improvement1,noise1,repeat,KF_type)
%% Problem setup
delta_t=0.01;
len=1;
mass=1;
g=9.8;
Tlimit=1;
sig_F=1e-3;
Q=[(sig_F*delta_t/mass/len^2)^2 0;0 0];
sitadot_rmse=0;
sita_rmse=0;
rng(123)
Runtime=0;
for iii=1:repeat
    states_true=[0;pi/4];
    sig_init=[pi/18;pi/18];
    Variance=[(sig_init(1))^2 0; 0 (sig_init(2))^2];
    states_est=states_true+normrnd(0, sig_init);
    %% generate profile
    Measurementnoise=noise1^2;
    t=delta_t:delta_t:Tlimit;
    v_true=zeros(1,length(t));
    v=v_true+(sig_F*delta_t/mass/len^2)*randn(1,length(t));
    stepnum=length(t);
    %v=v+randn(1,stepnum+1)*noise1;
    for i=1:stepnum
        states=states_true(:,i);
        states=states+[-g/len*sin(states(2))+v_true(i)/mass/len^2;states(1)]*delta_t;
        states_true=[states_true states];
    end
    sitadot=states_true(1,:);
    sita=states_true(2,:);
    Measurements=(mass*g*cos(sita)+mass*len*sitadot.^2).*sin(sita);
    Measurements_noisy=Measurements+randn(1,stepnum+1)*noise1;
    %% KF
    tic;
    for i=1:stepnum
        if(KF_type==0)
            states=states_est(:,i);
            F0=[1 -g/len*cos(states(2))*delta_t;delta_t 1];
            states=states+[-g/len*sin(states(2))+v(i)/mass/len^2;states(1)]*delta_t;
            Variance=F0*Variance*F0.'+Q;
            residual=Measurements_noisy(:,i+1)-(mass*g*cos(states(2))+mass*len*states(1)^2)*sin(states(2));
            H=[2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
            S=H*Variance*H.'+Measurementnoise;
            K=Variance*H.'*S^(-1);
            states0=states;
            states=states0+K*residual;
            if(improvement1==1)
                H2 = [2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
                temp=Variance-K*H2*Variance-Variance*H2'*K'+K*(H2*Variance*H2'+Measurementnoise)*K';       
                if(sum(diag(temp))<sum(diag(Variance)))
                    Variance=temp;
                else
                    states=states0;
                end
            else
                Variance=Variance-K*H*Variance;
            end
            states_est(:,i+1)=states;
        elseif(KF_type==0.5)
            states=states_est(:,i);
            F0=[1 -g/len*cos(states(2))*delta_t;delta_t 1];
            states=states+[-g/len*sin(states(2))+v(i)/mass/len^2;states(1)]*delta_t;
            Variance=F0*Variance*F0.'+Q;
            residual=Measurements_noisy(:,i+1)-(mass*g*cos(states(2))+mass*len*states(1)^2)*sin(states(2));
            H=[2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
            S=H*Variance*H.'+Measurementnoise;
            K=Variance*H.'*S^(-1);
            states0=states;
            states=states0+K*residual;
            steplennow = norm(K*residual);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                residual = Measurements_noisy(:,i+1)-(mass*g*cos(states(2))+mass*len*states(1)^2)*sin(states(2));
                H2 = [2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
                S2 = H2*Variance*H2.'+Measurementnoise;
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
                H2 = [2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
                temp=Variance-K*H2*Variance-Variance*H2'*K'+K*(H2*Variance*H2'+Measurementnoise)*K';       
                if(sum(diag(temp))<sum(diag(Variance)))
                    Variance=temp;
                else
                    states=states0;
                end
            else
                Variance=Variance-K*H*Variance;
            end
            states_est(:,i+1)=states;
        elseif(KF_type==1)
            L=real(sqrtm(Variance));
            state=states_est(:,i);
            n=2;
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
                states(:,ii)=states(:,ii)+[-g/len*sin(states(2,ii))+v(i)/mass/len^2;states(1,ii)]*delta_t;
                if(ii==1)
                    state=states(:,ii)*lambda/(n+lambda);
                else
                    state=state+states(:,ii)/(n+lambda)/2;
                end
            end
            Variance=(states(:,1)-state)*(states(:,1)-state).'*(lambda/(n+lambda)+1-alpha^2+beta)+Q;
            for ii=2:2*n+1
                Variance=Variance+(state-states(:,ii))*(state-states(:,ii)).'/(n+lambda)/2;
            end
            % obtain measurement
            z = Measurements_noisy(:,i+1);
            % Predict Measurement From Propagated Sigma Points
            measures=zeros(mnum,2*n+1);
            for ii=1:2*n+1
                measures(:,ii) = (mass*g*cos(states(2,ii))+mass*len*states(1,ii)^2)*sin(states(2,ii)); % predicted measurement
                if(ii==1)
                    m_exp=lambda/(n+lambda)*measures(:,ii);
                else
                    m_exp=m_exp+1/(n+lambda)/2*measures(:,ii);
                end
            end
            % Estimate Mean And Covariance of Predicted Measurement
            Py=(lambda/(n+lambda)+1-alpha^2+beta)*(measures(:,1)-m_exp)*(measures(:,1)-m_exp).'+Measurementnoise;
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
                    measures(:,ii)=(mass*g*cos(states(2,ii))+mass*len*states(1,ii)^2)*sin(states(2,ii));
                    if(ii==1)
                        m_exp=lambda/(n+lambda)*measures(:,ii);
                    else
                        m_exp=m_exp+1/(n+lambda)/2*measures(:,ii);
                    end
                end
                % Estimate Mean And Covariance of Predicted Measurement
                Py=(lambda/(n+lambda)+1-alpha^2+beta)*(measures(:,1)-m_exp)*(measures(:,1)-m_exp).'+Measurementnoise;
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
        elseif(KF_type==2)
            %CKF
            L=real(sqrtm(Variance));
            state=states_est(:,i);
            n=2;
            mnum=1;
            states=zeros(n,n*2);
            for ii=1:n
                states(:,ii)=state+sqrt(n)*L(:,ii);
                states(:,n+ii)=state-sqrt(n)*L(:,ii);  
            end
            state=0;
            for ii=1:2*n
                states(:,ii)=states(:,ii)+[-g/len*sin(states(2,ii))+v(i)/mass/len^2;states(1,ii)]*delta_t;
                state=state+states(:,ii)/n/2;
            end
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
                measures(:,ii) = (mass*g*cos(states(2,ii))+mass*len*states(1,ii)^2)*sin(states(2,ii)); % predicted measurement
                m_exp=m_exp+1/n/2*measures(:,ii);
            end
            % Estimate Mean And Covariance of Predicted Measurement
            Py=Measurementnoise;
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
                    measures(:,ii) = (mass*g*cos(states(2,ii))+mass*len*states(1,ii)^2)*sin(states(2,ii)); % predicted measurement
                    m_exp=m_exp+1/n/2*measures(:,ii);
                end
                % Estimate Mean And Covariance of Predicted Measurement
                Py=Measurementnoise;
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
        elseif(KF_type==3)
            %2-EKF
            states=states_est(:,i);
            n=length(states);
            mnum=1;
            F=[1 -g/len*cos(states(2))*delta_t;delta_t 1];
            states=states+[-g/len*sin(states(2))+v(i)/mass/len^2;states(1)]*delta_t;
            Fxx=zeros(n,n,n);
            deltaP=zeros(n,n);
            Fxx(2,2,1)=g/len*sin(states(2))*delta_t;
            for ii=1:n
                states(ii)=states(ii)+0.5*sum(diag(Fxx(:,:,ii)*Variance));
                for jj=1:n
                    deltaP(ii,jj)=sum(diag(Fxx(:,:,ii)*Variance*Fxx(:,:,jj)*Variance));
                end
            end
            Variance=F*Variance*F.'+0.5*deltaP+Q;
            %update
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=[2*mass*len*sin(states(2)) 2*mass*len*states(1)*cos(states(2));
                2*mass*len*states(1)*cos(states(2)) -2*mass*g*sin(2*states(2))-mass*len*states(1)^2*sin(states(2))]; 
            H=[2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
            S=H*Variance*H.'+Measurementnoise;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*Variance*Hxx(:,:,jj)*Variance));
                end
            end
            K=Variance*H.'*S^(-1);
            states0=states;
            residual=Measurements_noisy(:,i+1)-(mass*g*cos(states(2))+mass*len*states(1)^2)*sin(states(2));
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));
            end
            states=states0+K*residual;
            if(improvement1==1)
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[2*mass*len*sin(states(2)) 2*mass*len*states(1)*cos(states(2));
                2*mass*len*states(1)*cos(states(2)) -2*mass*g*sin(2*states(2))-mass*len*states(1)^2*sin(states(2))]; 
                H2=[2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
                S2=H2*Variance*H2.'+Measurementnoise;
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
            n=length(states);
            mnum=1;
            F=[1 -g/len*cos(states(2))*delta_t;delta_t 1];
            states=states+[-g/len*sin(states(2))+v(i)/mass/len^2;states(1)]*delta_t;
            Fxx=zeros(n,n,n);
            deltaP=zeros(n,n);
            Fxx(2,2,1)=g/len*sin(states(2))*delta_t;
            for ii=1:n
                states(ii)=states(ii)+0.5*sum(diag(Fxx(:,:,ii)*Variance));
                for jj=1:n
                    deltaP(ii,jj)=sum(diag(Fxx(:,:,ii)*Variance*Fxx(:,:,jj)*Variance));
                end
            end
            Variance=F*Variance*F.'+0.5*deltaP+Q;
            %update
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=[2*mass*len*sin(states(2)) 2*mass*len*states(1)*cos(states(2));
                2*mass*len*states(1)*cos(states(2)) -2*mass*g*sin(2*states(2))-mass*len*states(1)^2*sin(states(2))]; 
            H=[2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
            S=H*Variance*H.'+Measurementnoise;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*Variance*Hxx(:,:,jj)*Variance));
                end
            end
            K=Variance*H.'*S^(-1);
            states0=states;
            residual=Measurements_noisy(:,i+1)-(mass*g*cos(states(2))+mass*len*states(1)^2)*sin(states(2));
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));
            end
            states=states0+K*residual;
            steplennow = norm(K*residual);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[2*mass*len*sin(states(2)) 2*mass*len*states(1)*cos(states(2));
                2*mass*len*states(1)*cos(states(2)) -2*mass*g*sin(2*states(2))-mass*len*states(1)^2*sin(states(2))]; 
                H2 = [2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
                residual = Measurements_noisy(:,i+1)-(mass*g*cos(states(2))+mass*len*states(1)^2)*sin(states(2));
                for ii=1:mnum
                    residual(ii)=residual(ii)-0.5*(states0-states).'*Hxx2(:,:,ii)*(states0-states);
                    residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));
                end
                S2 = H2*Variance*H2.'+Measurementnoise;
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
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[2*mass*len*sin(states(2)) 2*mass*len*states(1)*cos(states(2));
                2*mass*len*states(1)*cos(states(2)) -2*mass*g*sin(2*states(2))-mass*len*states(1)^2*sin(states(2))]; 
                H2=[2*mass*len*states(1)*sin(states(2)) mass*g*cos(2*states(2))+mass*len*states(1)^2*cos(states(2))];
                S2=H2*Variance*H2.'+Measurementnoise;
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
    sitadot_est=states_est(1,:);
    sita_est=states_est(2,:);
    sitadot_rmse=sitadot_rmse+(sitadot(length(sitadot))-sitadot_est(length(sitadot)))^2/repeat;
    sita_rmse=sita_rmse+(sita(length(sitadot))-sita_est(length(sitadot)))^2/repeat;
end
sita_rmse=sqrt(sita_rmse);
sitadot_rmse=sqrt(sitadot_rmse);
end