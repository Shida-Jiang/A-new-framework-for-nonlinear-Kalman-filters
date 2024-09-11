% Application 3: Synchronous Generator State Estimation
close all
clc
clear
repeat=10000; %The number of times that the simulation is repeated at each measurement noise setup
%This simulation takes about 30 minutes, you can reduce the repeat times to reduce the run time
%%
noise1_ref=1;
noise2_ref=1;
scale=-6:0.5:-1; % The range of measurement noise is 10^scale
x1plotEKF_1=[];
x2plotEKF_1=[];
x3plotEKF_1=[];
x4plotEKF_1=[];
Time_EKF_1=[];
x1plotEKF2_1=[];
x2plotEKF2_1=[];
x3plotEKF2_1=[];
x4plotEKF2_1=[];
Time_EKF2_1=[];
x1plotUKF_1=[];
x2plotUKF_1=[];
x3plotUKF_1=[];
x4plotUKF_1=[];
Time_UKF_1=[];
x1plotCKF_1=[];
x2plotCKF_1=[];
x3plotCKF_1=[];
x4plotCKF_1=[];
Time_CKF_1=[];

x1plotEKF_2=[];
x2plotEKF_2=[];
x3plotEKF_2=[];
x4plotEKF_2=[];
Time_EKF_2=[];
x1plotEKF2_2=[];
x2plotEKF2_2=[];
x3plotEKF2_2=[];
x4plotEKF2_2=[];
Time_EKF2_2=[];
x1plotUKF_2=[];
x2plotUKF_2=[];
x3plotUKF_2=[];
x4plotUKF_2=[];
Time_UKF_2=[];
x1plotCKF_2=[];
x2plotCKF_2=[];
x3plotCKF_2=[];
x4plotCKF_2=[];
Time_CKF_2=[];
x1plotIEKF_1=[];
x2plotIEKF_1=[];
x3plotIEKF_1=[];
x4plotIEKF_1=[];
Time_IEKF_1=[];
x1plotIEKF_2=[];
x2plotIEKF_2=[];
x3plotIEKF_2=[];
x4plotIEKF_2=[];
Time_IEKF_2=[];
x1plotIEKF2_1=[];
x2plotIEKF2_1=[];
x3plotIEKF2_1=[];
x4plotIEKF2_1=[];
Time_IEKF2_1=[];
x1plotIEKF2_2=[];
x2plotIEKF2_2=[];
x3plotIEKF2_2=[];
x4plotIEKF2_2=[];
Time_IEKF2_2=[];
for i=scale
    noise1=noise1_ref*10^i;
    noise2=noise2_ref*10^i;

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(0,noise1,repeat,0.5);
    x1plotIEKF_1=[x1plotIEKF_1 x1_rmse];
    x2plotIEKF_1=[x2plotIEKF_1 x2_rmse];
    x3plotIEKF_1=[x3plotIEKF_1 x3_rmse];
    x4plotIEKF_1=[x4plotIEKF_1 x4_rmse];
    Time_IEKF_1=[Time_IEKF_1 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(1,noise1,repeat,0.5);
    x1plotIEKF_2=[x1plotIEKF_2 x1_rmse];
    x2plotIEKF_2=[x2plotIEKF_2 x2_rmse];
    x3plotIEKF_2=[x3plotIEKF_2 x3_rmse];
    x4plotIEKF_2=[x4plotIEKF_2 x4_rmse];
    Time_IEKF_2=[Time_IEKF_2 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(1,noise1,repeat,0);
    x1plotEKF_2=[x1plotEKF_2 x1_rmse];
    x2plotEKF_2=[x2plotEKF_2 x2_rmse];
    x3plotEKF_2=[x3plotEKF_2 x3_rmse];
    x4plotEKF_2=[x4plotEKF_2 x4_rmse];
    Time_EKF_2=[Time_EKF_2 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(0,noise1,repeat,0);
    x1plotEKF_1=[x1plotEKF_1 x1_rmse];
    x2plotEKF_1=[x2plotEKF_1 x2_rmse];
    x3plotEKF_1=[x3plotEKF_1 x3_rmse];
    x4plotEKF_1=[x4plotEKF_1 x4_rmse];
    Time_EKF_1=[Time_EKF_1 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(0,noise1,repeat,3);
    x1plotEKF2_1=[x1plotEKF2_1 x1_rmse];
    x2plotEKF2_1=[x2plotEKF2_1 x2_rmse];
    x3plotEKF2_1=[x3plotEKF2_1 x3_rmse];
    x4plotEKF2_1=[x4plotEKF2_1 x4_rmse];
    Time_EKF2_1=[Time_EKF2_1 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(1,noise1,repeat,3);
    x1plotEKF2_2=[x1plotEKF2_2 x1_rmse];
    x2plotEKF2_2=[x2plotEKF2_2 x2_rmse];
    x3plotEKF2_2=[x3plotEKF2_2 x3_rmse];
    x4plotEKF2_2=[x4plotEKF2_2 x4_rmse];
    Time_EKF2_2=[Time_EKF2_2 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(0,noise1,repeat,3.5);
    x1plotIEKF2_1=[x1plotIEKF2_1 x1_rmse];
    x2plotIEKF2_1=[x2plotIEKF2_1 x2_rmse];
    x3plotIEKF2_1=[x3plotIEKF2_1 x3_rmse];
    x4plotIEKF2_1=[x4plotIEKF2_1 x4_rmse];
    Time_IEKF2_1=[Time_IEKF2_1 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(1,noise1,repeat,3.5);
    x1plotIEKF2_2=[x1plotIEKF2_2 x1_rmse];
    x2plotIEKF2_2=[x2plotIEKF2_2 x2_rmse];
    x3plotIEKF2_2=[x3plotIEKF2_2 x3_rmse];
    x4plotIEKF2_2=[x4plotIEKF2_2 x4_rmse];
    Time_IEKF2_2=[Time_IEKF2_2 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(0,noise1,repeat,1);
    x1plotUKF_1=[x1plotUKF_1 x1_rmse];
    x2plotUKF_1=[x2plotUKF_1 x2_rmse];
    x3plotUKF_1=[x3plotUKF_1 x3_rmse];
    x4plotUKF_1=[x4plotUKF_1 x4_rmse];
    Time_UKF_1=[Time_UKF_1 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(1,noise1,repeat,1);
    x1plotUKF_2=[x1plotUKF_2 x1_rmse];
    x2plotUKF_2=[x2plotUKF_2 x2_rmse];
    x3plotUKF_2=[x3plotUKF_2 x3_rmse];
    x4plotUKF_2=[x4plotUKF_2 x4_rmse];
    Time_UKF_2=[Time_UKF_2 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(0,noise1,repeat,2);
    x1plotCKF_1=[x1plotCKF_1 x1_rmse];
    x2plotCKF_1=[x2plotCKF_1 x2_rmse];
    x3plotCKF_1=[x3plotCKF_1 x3_rmse];
    x4plotCKF_1=[x4plotCKF_1 x4_rmse];
    Time_CKF_1=[Time_CKF_1 Runtime];

    [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(1,noise1,repeat,2);
    x1plotCKF_2=[x1plotCKF_2 x1_rmse];
    x2plotCKF_2=[x2plotCKF_2 x2_rmse];
    x3plotCKF_2=[x3plotCKF_2 x3_rmse];
    x4plotCKF_2=[x4plotCKF_2 x4_rmse];
    Time_CKF_2=[Time_CKF_2 Runtime];
    
end

%% display Runtime
disp(['EKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_EKF_1)), ' ms'])
disp(['EKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_EKF_2)), ' ms'])
disp(['IEKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_IEKF_1)), ' ms'])
% disp(['IEKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_IEKF_2)), ' ms'])
disp(['EKF2 average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_EKF2_1)), ' ms'])
disp(['EKF2 average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_EKF2_2)), ' ms'])
% disp(['IEKF2 average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_IEKF2_1)), ' ms'])
% disp(['IEKF2 average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_IEKF2_2)), ' ms'])
disp(['UKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_UKF_1)), ' ms'])
disp(['UKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_UKF_2)), ' ms'])
disp(['CKF average runtime (old framework): ',sprintf('%.2f', 1000*mean(Time_CKF_1)), ' ms'])
disp(['CKF average runtime (new framework): ',sprintf('%.2f', 1000*mean(Time_CKF_2)), ' ms'])
%% figures
f1=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
h1_2=plot(10.^scale, x1plotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, x1plotEKF_1, '--o', 'DisplayName', "EKF (old)", 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, x1plotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, x1plotEKF2_1, '--s', 'DisplayName', "EKF2 (old)", 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, x1plotUKF_2, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, x1plotUKF_1, '--^', 'DisplayName', "UKF (old)", 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, x1plotCKF_2, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, x1plotCKF_1, '--x', 'DisplayName', "CKF (old)", 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, x1plotIEKF_2, '-d', 'Color', '#4DBEEE','LineWidth',1);
h5_1=plot(10.^scale, x1plotIEKF_1, '--d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
%h6_2=plot(10.^scale, x1plotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, x1plotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
ylabel('Angle RMSE (rad)',FontSize=12)
set(gca, 'YTick', [0.0001 0.001 0.01]);
ylim([1e-4 1e-2])
xlim([10^-6 0.1])
h1_1=plot(nan, nan, '--', 'Color', 'k', 'DisplayName', 'Old framework','LineWidth',1);
h2_1=plot(nan, nan, 'Color', "k", 'DisplayName', 'New framework','LineWidth',1);
leg=legend([h1_1 h2_1], 'Location','southeast',FontSize=11);
title(leg,'Line styles')
grid on

nexttile
hold on
h1_2=plot(10.^scale, x2plotEKF_2, '-o', 'DisplayName', "EKF (new)", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, x2plotEKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, x2plotEKF2_2, '-s', 'DisplayName', "EKF2 (new)", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, x2plotEKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, x2plotUKF_2, '-^', 'DisplayName', "UKF (new)", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, x2plotUKF_1, '--^', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, x2plotCKF_2, '-x', 'DisplayName', "CKF (new)", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, x2plotCKF_1, '--x', 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, x2plotIEKF_2, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
h5_1=plot(10.^scale, x2plotIEKF_1, '--d', 'Color', '#4DBEEE','LineWidth',1);
%h6_2=plot(10.^scale, x2plotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, x2plotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
h1_1=plot(nan, nan, 'o', 'Color', '#0072BD', 'DisplayName', 'EKF','LineWidth',1);
h2_1=plot(nan, nan, 's', 'Color', "#D95319", 'DisplayName', 'EKF2','LineWidth',1);
h3_1=plot(nan, nan, '^', 'Color', "#EDB120", 'DisplayName', 'UKF','LineWidth',1);
h4_1=plot(nan, nan, 'x', 'Color', "#7E2F8E", 'DisplayName', 'CKF','LineWidth',1);
h5_1=plot(nan, nan, 'd', 'Color', "#4DBEEE", 'DisplayName', 'IEKF','LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
xlabel('Measurement standard deviations (per unit)',FontSize=12)
ylabel('Rotor speed RMSE (rad/s)',FontSize=12)
ylim([1e-2 1e-1])
xlim([10^-6 0.1])
%set(gca, 'YTick', [0.001 0.01 0.1 1]);
leg=legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast','NumColumns',1,FontSize=11);
title(leg,'Line colors')
grid on
%exportgraphics(f1,'Synchronous-Generator.png','Resolution',900)

% %%
% figure
% f2=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
% nexttile
% hold on
% h1_2=plot(10.^scale, x1plotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, x1plotEKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(10.^scale, x1plotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, x1plotEKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(10.^scale, x1plotIEKF_2, '-d', 'Color', '#4DBEEE','LineWidth',1);
% h5_1=plot(10.^scale, x1plotIEKF_1, '--d', 'Color', '#4DBEEE','LineWidth',1);
% h6_2=plot(10.^scale, x1plotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, x1plotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
% ylabel('Angle RMSE (rad)',FontSize=12)
% set(gca, 'YTick', [0.0001 0.001 0.01]);
% ylim([1e-4 1e-2])
% xlim([10^-6 0.1])
% grid on
% %legend([h1_2 h2_2 h3_2 h4_2], 'Location','southeast',FontSize=10)
% 
% % Summary
% % disp(append("X-axis RMS error =",num2str(x_rmse)))
% % disp(append("Y-axis RMS error =",num2str(y_rmse)))
% % disp(append("Angle RMS error =",num2str(phi_rmse)," degree"))
% 
% nexttile
% hold on
% h1_2=plot(10.^scale, x2plotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, x2plotEKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(10.^scale, x2plotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, x2plotEKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(10.^scale, x2plotIEKF_2, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
% h5_1=plot(10.^scale, x2plotIEKF_1, '--d', 'Color', '#4DBEEE','LineWidth',1);
% h6_2=plot(10.^scale, x2plotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, x2plotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
% xlabel('Measurement standard deviations (per unit)',FontSize=12)
% ylabel('Rotor speed RMSE (rad/s)',FontSize=12)
% ylim([1e-2 1e-1])
% xlim([10^-6 0.1])
% %set(gca, 'YTick', [0.001 0.01 0.1 1]);
% legend([h1_2 h2_2 h5_2 h6_2], 'Location','southeast',FontSize=12)
% grid on
% exportgraphics(f2,'Synchronous-Generator-IKF.png','Resolution',900)
function [x1_rmse,x2_rmse,x3_rmse,x4_rmse,Runtime]=simulation(improvement1,noise1,repeat,KFtype)
%% Problem setup
delta_t=1e-4;
w0=2*60*pi;
J=13;
x_d_dot=0.375;
x_d=2.06;
x_q_dot=0.375;
x_q=1.214;
T_do=0.131;
T_qo=0.0131;
D=0.05;
Tlimit=1e-2;
x1_rmse=0;
x2_rmse=0;
x3_rmse=0;
x4_rmse=0;
process_sigma=[1e-5;1e-8;1e-5;1e-5];
Q=diag(process_sigma.^2);
sig_init = [0.01; 1e-5; 0.01; 0.01]; % standard deviation of initial guess
rng(123)
Runtime=0;
for iii=1:repeat
    states_true=[0.4;0;0;0];
    states_est=states_true+normrnd(0, sig_init);
    Variance=diag(sig_init.^2);
    %% generate profile
    Measurementnoise=noise1^2;
    t=delta_t:delta_t:Tlimit;
    v=zeros(3,length(t));
    v(1,:)=0.8*ones(1,length(t));
    v(2,:)=2.11+2*(t-delta_t);
    v(3,:)=1.002*ones(1,length(t));
    stepnum=length(t);
    Measurements=[];
    for i=1:stepnum
        states=states_true(:,i);
        states=states+delta_t*[w0*states(2);
            (v(1,i)-D*states(2)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*sin(2*states(1))*(1/x_q-1/x_d_dot)))/J;
            (v(2,i)-states(3)-(x_d-x_d_dot)*(states(3)-v(3,i)*cos(states(1)))/x_d_dot)/T_do;
            (-states(4)+(x_q-x_q_dot)*v(3,i)*sin(states(1))/x_q)/T_qo]+process_sigma.*randn(4,1);
        states_true=[states_true states];
        measurement=v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1));
        Measurements=[Measurements measurement];
    end
    x=states_true(1,:);
    y=states_true(2,:);
    Measurements = [0 Measurements];
    Measurements_noisy=Measurements+noise1*randn(1,stepnum+1);
    %% KF
    tic;
    for i=1:stepnum
        if(KFtype==0)
            states=states_est(:,i);
            F=diag([1 1 1 1])+delta_t*[0 w0 0 0;
                -(v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*cos(2*states(1))*(1/x_q-1/x_d_dot))/J -D/J -v(3,i)/x_d_dot*sin(states(1))/J 0;
                -(x_d-x_d_dot)*v(3,i)*sin(states(1))/x_d_dot/T_do 0 -(1+(x_d-x_d_dot)/x_d_dot)/T_do 0;
                (x_q-x_q_dot)*v(3,i)*cos(states(1))/x_q/T_qo 0 0 -1/T_qo];
            states=states+delta_t*[w0*states(2);
            (v(1,i)-D*states(2)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*sin(2*states(1))*(1/x_q-1/x_d_dot)))/J;
            (v(2,i)-states(3)-(x_d-x_d_dot)*(states(3)-v(3,i)*cos(states(1)))/x_d_dot)/T_do;
            (-states(4)+(x_q-x_q_dot)*v(3,i)*sin(states(1))/x_q)/T_qo];
            Variance=F*Variance*F.'+Q;
            residual=Measurements_noisy(:,i+1)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1)));
            H=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
            S=H*Variance*H.'+Measurementnoise;
            K=Variance*H.'*S^(-1);
            states0=states;
            states=states+K*residual;
            if(improvement1==1)
                H2=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
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
        elseif(KFtype==0.5)
            states=states_est(:,i);
            F=diag([1 1 1 1])+delta_t*[0 w0 0 0;
                -(v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*cos(2*states(1))*(1/x_q-1/x_d_dot))/J -D/J -v(3,i)/x_d_dot*sin(states(1))/J 0;
                -(x_d-x_d_dot)*v(3,i)*sin(states(1))/x_d_dot/T_do 0 -(1+(x_d-x_d_dot)/x_d_dot)/T_do 0;
                (x_q-x_q_dot)*v(3,i)*cos(states(1))/x_q/T_qo 0 0 -1/T_qo];
            states=states+delta_t*[w0*states(2);
            (v(1,i)-D*states(2)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*sin(2*states(1))*(1/x_q-1/x_d_dot)))/J;
            (v(2,i)-states(3)-(x_d-x_d_dot)*(states(3)-v(3,i)*cos(states(1)))/x_d_dot)/T_do;
            (-states(4)+(x_q-x_q_dot)*v(3,i)*sin(states(1))/x_q)/T_qo];
            Variance=F*Variance*F.'+Q;
            residual=Measurements_noisy(:,i+1)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1)));
            H=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
            S=H*Variance*H.'+Measurementnoise;
            K=Variance*H.'*S^(-1);
            states0=states;
            states=states+K*residual;
            steplennow = norm(K*residual);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                residual = Measurements_noisy(:,i+1)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1)));
                H2 = [v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
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
                H2=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
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
        elseif(KFtype==1)
            L=real(sqrtm(Variance));
            state=states_est(:,i);
            n=4;
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
                states(:,ii)=states(:,ii)+delta_t*[w0*states(2,ii);
                    (v(1,i)-D*states(2,ii)-(v(3,i)/x_d_dot*states(3,ii)*sin(states(1,ii))+v(3,i)^2/2*sin(2*states(1,ii))*(1/x_q-1/x_d_dot)))/J;
                    (v(2,i)-states(3,ii)-(x_d-x_d_dot)*(states(3,ii)-v(3,i)*cos(states(1,ii)))/x_d_dot)/T_do;
                    (-states(4,ii)+(x_q-x_q_dot)*v(3,i)*sin(states(1,ii))/x_q)/T_qo];
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
                measures(:,ii) = (v(3,i)/x_d_dot*states(3,ii)*sin(states(1,ii))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1,ii))); % predicted measurement
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
            state0=state;
            state=state+K*(z-m_exp);
            %update
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
                    measures(:,ii) = (v(3,i)/x_d_dot*states(3,ii)*sin(states(1,ii))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1,ii))); % predicted measurement
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
        elseif(KFtype==2)
            %CKF
            L=real(sqrtm(Variance));
            state=states_est(:,i);
            n=4;
            mnum=1;
            states=zeros(n,n*2);
            for ii=1:n
                states(:,ii)=state+sqrt(n)*L(:,ii);
                states(:,n+ii)=state-sqrt(n)*L(:,ii);  
            end
            state=0;
            for ii=1:2*n
                states(:,ii)=states(:,ii)+delta_t*[w0*states(2,ii);
                    (v(1,i)-D*states(2,ii)-(v(3,i)/x_d_dot*states(3,ii)*sin(states(1,ii))+v(3,i)^2/2*sin(2*states(1,ii))*(1/x_q-1/x_d_dot)))/J;
                    (v(2,i)-states(3,ii)-(x_d-x_d_dot)*(states(3,ii)-v(3,i)*cos(states(1,ii)))/x_d_dot)/T_do;
                    (-states(4,ii)+(x_q-x_q_dot)*v(3,i)*sin(states(1,ii))/x_q)/T_qo];
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
                measures(:,ii) = (v(3,i)/x_d_dot*states(3,ii)*sin(states(1,ii))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1,ii))); % predicted measurement
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
                    measures(:,ii) = (v(3,i)/x_d_dot*states(3,ii)*sin(states(1,ii))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1,ii))); % predicted measurement
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
            
        elseif(KFtype==3)
            %2-EKF
            states=states_est(:,i);
            states=states+delta_t*[w0*states(2);
            (v(1,i)-D*states(2)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*sin(2*states(1))*(1/x_q-1/x_d_dot)))/J;
            (v(2,i)-states(3)-(x_d-x_d_dot)*(states(3)-v(3,i)*cos(states(1)))/x_d_dot)/T_do;
            (-states(4)+(x_q-x_q_dot)*v(3,i)*sin(states(1))/x_q)/T_qo];
            n=length(states);
            Fxx=zeros(n,n,n);
            deltaP=zeros(n,n);
            Fxx(1,1,2)=(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2*2*sin(2*states(1))*(1/x_q-1/x_d_dot))/J*delta_t;
            Fxx(1,3,2)=-v(3,i)/x_d_dot*cos(states(1))/J*delta_t;
            Fxx(3,1,2)=-v(3,i)/x_d_dot*cos(states(1))/J*delta_t;
            Fxx(1,1,3)=-(x_d-x_d_dot)*v(3,i)*cos(states(1))/x_d_dot/T_do*delta_t;
            Fxx(1,1,4)=-(x_q-x_q_dot)*v(3,i)*sin(states(1))/x_q/T_qo*delta_t;
            for ii=1:n
                states(ii)=states(ii)+0.5*sum(diag(Fxx(:,:,ii)*Variance));
                for jj=1:n
                    deltaP(ii,jj)=sum(diag(Fxx(:,:,ii)*Variance*Fxx(:,:,jj)*Variance));
                end
            end
            F=diag([1 1 1 1])+delta_t*[0 w0 0 0;
                -(v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*cos(2*states(1))*(1/x_q-1/x_d_dot))/J -D/J -v(3,i)/x_d_dot*sin(states(1))/J 0;
                -(x_d-x_d_dot)*v(3,i)*sin(states(1))/x_d_dot/T_do 0 -(1+(x_d-x_d_dot)/x_d_dot)/T_do 0;
                (x_q-x_q_dot)*v(3,i)*cos(states(1))/x_q/T_qo 0 0 -1/T_qo];
            Variance=F*Variance*F.'+0.5*deltaP+Q;
            %update
            mnum=1;
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=[-v(3,i)/x_d_dot*states(3)*sin(states(1))-2*v(3,i)^2*(1/x_q-1/x_d_dot)*sin(2*states(1)) 0 v(3,i)/x_d_dot*cos(states(1)) 0;
                0 0 0 0;
                v(3,i)/x_d_dot*cos(states(1)) 0 0 0;
                0 0 0 0];
            H=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
            S=H*Variance*H.'+Measurementnoise;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*Variance*Hxx(:,:,jj)*Variance));
                end
            end
            K=Variance*H.'*S^(-1);
            states0=states;
            residual=Measurements_noisy(:,i+1)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1)));
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));
            end
            states=states+K*residual;
            if(improvement1==1)
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[-v(3,i)/x_d_dot*states(3)*sin(states(1))-2*v(3,i)^2*(1/x_q-1/x_d_dot)*sin(2*states(1)) 0 v(3,i)/x_d_dot*cos(states(1)) 0;
                    0 0 0 0;
                    v(3,i)/x_d_dot*cos(states(1)) 0 0 0;
                    0 0 0 0];
                H2=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
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
            states=states+delta_t*[w0*states(2);
            (v(1,i)-D*states(2)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*sin(2*states(1))*(1/x_q-1/x_d_dot)))/J;
            (v(2,i)-states(3)-(x_d-x_d_dot)*(states(3)-v(3,i)*cos(states(1)))/x_d_dot)/T_do;
            (-states(4)+(x_q-x_q_dot)*v(3,i)*sin(states(1))/x_q)/T_qo];
            n=length(states);
            Fxx=zeros(n,n,n);
            deltaP=zeros(n,n);
            Fxx(1,1,2)=(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2*2*sin(2*states(1))*(1/x_q-1/x_d_dot))/J*delta_t;
            Fxx(1,3,2)=-v(3,i)/x_d_dot*cos(states(1))/J*delta_t;
            Fxx(3,1,2)=-v(3,i)/x_d_dot*cos(states(1))/J*delta_t;
            Fxx(1,1,3)=-(x_d-x_d_dot)*v(3,i)*cos(states(1))/x_d_dot/T_do*delta_t;
            Fxx(1,1,4)=-(x_q-x_q_dot)*v(3,i)*sin(states(1))/x_q/T_qo*delta_t;
            for ii=1:n
                states(ii)=states(ii)+0.5*sum(diag(Fxx(:,:,ii)*Variance));
                for jj=1:n
                    deltaP(ii,jj)=sum(diag(Fxx(:,:,ii)*Variance*Fxx(:,:,jj)*Variance));
                end
            end
            F=diag([1 1 1 1])+delta_t*[0 w0 0 0;
                -(v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*cos(2*states(1))*(1/x_q-1/x_d_dot))/J -D/J -v(3,i)/x_d_dot*sin(states(1))/J 0;
                -(x_d-x_d_dot)*v(3,i)*sin(states(1))/x_d_dot/T_do 0 -(1+(x_d-x_d_dot)/x_d_dot)/T_do 0;
                (x_q-x_q_dot)*v(3,i)*cos(states(1))/x_q/T_qo 0 0 -1/T_qo];
            Variance=F*Variance*F.'+0.5*deltaP+Q;
            %update
            mnum=1;
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=[-v(3,i)/x_d_dot*states(3)*sin(states(1))-2*v(3,i)^2*(1/x_q-1/x_d_dot)*sin(2*states(1)) 0 v(3,i)/x_d_dot*cos(states(1)) 0;
                0 0 0 0;
                v(3,i)/x_d_dot*cos(states(1)) 0 0 0;
                0 0 0 0];
            H=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
            S=H*Variance*H.'+Measurementnoise;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*Variance*Hxx(:,:,jj)*Variance));
                end
            end
            K=Variance*H.'*S^(-1);
            states0=states;
            residual=Measurements_noisy(:,i+1)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1)));
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));
            end
            states=states+K*residual;
            steplennow = norm(K*residual);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                residual = Measurements_noisy(:,i+1)-(v(3,i)/x_d_dot*states(3)*sin(states(1))+v(3,i)^2/2*(1/x_q-1/x_d_dot)*sin(2*states(1)));
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=[-v(3,i)/x_d_dot*states(3)*sin(states(1))-2*v(3,i)^2*(1/x_q-1/x_d_dot)*sin(2*states(1)) 0 v(3,i)/x_d_dot*cos(states(1)) 0;
                    0 0 0 0;
                    v(3,i)/x_d_dot*cos(states(1)) 0 0 0;
                    0 0 0 0];
                for ii=1:mnum
                    residual(ii)=residual(ii)-0.5*(states0-states).'*Hxx2(:,:,ii)*(states0-states);
                    residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*Variance));  
                end
                H2 = [v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
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
                Hxx2(:,:,1)=[-v(3,i)/x_d_dot*states(3)*sin(states(1))-2*v(3,i)^2*(1/x_q-1/x_d_dot)*sin(2*states(1)) 0 v(3,i)/x_d_dot*cos(states(1)) 0;
                    0 0 0 0;
                    v(3,i)/x_d_dot*cos(states(1)) 0 0 0;
                    0 0 0 0];
                H2=[v(3,i)/x_d_dot*states(3)*cos(states(1))+v(3,i)^2*(1/x_q-1/x_d_dot)*cos(2*states(1)) 0 v(3,i)/x_d_dot*sin(states(1)) 0];
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
    x_est=states_est(1,:);
    y_est=states_est(2,:);
    x1_rmse=x1_rmse+(x(length(x))-x_est(length(x)))^2/repeat;
    x2_rmse=x2_rmse+(y(length(x))-y_est(length(x)))^2/repeat;
    x3_rmse=x3_rmse+(states_true(3,length(x))-states_est(3,length(x)))^2/repeat; 
end
x1_rmse=sqrt(x1_rmse);
x2_rmse=sqrt(x1_rmse);
x3_rmse=sqrt(x1_rmse);
end
