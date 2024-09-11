%%% This code is adapted from:
% Y. Kim and H. Bang, Introduction to Kalman Filter and Its Applications, InTechOpen, 2018

% Application 1: 3D target tracking
close all
clc
clear

repeat = 10000; %The number of times that the simulation is repeated at each measurement noise setup
%This simulation takes about 30 minutes, you can reduce the repeat times to reduce the run time

%% simulation
scale=-4:0.5:2; % The range of measurement noise is 10^scale
i0 = -2; % show how the estimation convergences when the std of the measurement noise is 10^-2 
%Defining some parameters
xplotEKF_1=[];
xerrorplotEKF_1=[];
vxplotEKF_1=[];
Time_EKF_1=[];
xplotEKF2_1=[];
xerrorplotEKF2_1=[];
vxplotEKF2_1=[];
Time_EKF2_1=[];
xplotUKF_1=[];
xerrorplotUKF_1=[];
vxplotUKF_1=[];
Time_UKF_1=[];
xplotCKF_1=[];
xerrorplotCKF_1=[];
vxplotCKF_1=[];
Time_CKF_1=[];

xplotEKF_2=[];
xerrorplotEKF_2=[];
vxplotEKF_2=[];
Time_EKF_2=[];
xplotEKF2_2=[];
xerrorplotEKF2_2=[];
vxplotEKF2_2=[];
Time_EKF2_2=[];
xplotUKF_2=[];
xerrorplotUKF_2=[];
vxplotUKF_2=[];
Time_UKF_2=[];
xplotCKF_2=[];
xerrorplotCKF_2=[];
vxplotCKF_2=[];
Time_CKF_2=[];
xplotIEKF_1=[];
xerrorplotIEKF_1=[];
vxplotIEKF_1=[];
Time_IEKF_1=[];
xplotIEKF_2=[];
xerrorplotIEKF_2=[];
vxplotIEKF_2=[];
Time_IEKF_2=[];
xplotIEKF2_1=[];
xerrorplotIEKF2_1=[];
vxplotIEKF2_1=[];
Time_IEKF2_1=[];
xplotIEKF2_2=[];
xerrorplotIEKF2_2=[];
vxplotIEKF2_2=[];
Time_IEKF2_2=[];
for i=scale
    magnitude=10^i;
    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(0,magnitude,repeat,0);
    len=length(x_RMSE(1,:));
    xplotEKF_1=[xplotEKF_1 x_RMSE(1,len)];
    xerrorplotEKF_1=[xerrorplotEKF_1 x_error_self_est];
    vxplotEKF_1=[vxplotEKF_1 x_RMSE(4,len)];
    Time_EKF_1=[Time_EKF_1 Runtime];
    if i==i0
        x_EKF_1=x_RMSE(1,:);
        y_EKF_1=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(1,magnitude,repeat,0);
    xplotEKF_2=[xplotEKF_2 x_RMSE(1,len)];
    xerrorplotEKF_2=[xerrorplotEKF_2 x_error_self_est];
    vxplotEKF_2=[vxplotEKF_2 x_RMSE(4,len)];
    Time_EKF_2=[Time_EKF_2 Runtime];
    if i==i0
        x_EKF_2=x_RMSE(1,:);
        y_EKF_2=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(0,magnitude,repeat,0.5);
    xplotIEKF_1=[xplotIEKF_1 x_RMSE(1,len)];
    xerrorplotIEKF_1=[xerrorplotIEKF_1 x_error_self_est];
    vxplotIEKF_1=[vxplotIEKF_1 x_RMSE(4,len)];
    Time_IEKF_1=[Time_IEKF_1 Runtime];
    if i==i0
        x_IEKF_1=x_RMSE(1,:);
        y_IEKF_1=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(1,magnitude,repeat,0.5);
    xplotIEKF_2=[xplotIEKF_2 x_RMSE(1,len)];
    xerrorplotIEKF_2=[xerrorplotIEKF_2 x_error_self_est];
    vxplotIEKF_2=[vxplotIEKF_2 x_RMSE(4,len)];
    Time_IEKF_2=[Time_IEKF_2 Runtime];
    if i==i0
        x_IEKF_2=x_RMSE(1,:);
        y_IEKF_2=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(0,magnitude,repeat,3);
    xplotEKF2_1=[xplotEKF2_1 x_RMSE(1,len)];
    xerrorplotEKF2_1=[xerrorplotEKF2_1 x_error_self_est];
    vxplotEKF2_1=[vxplotEKF2_1 x_RMSE(4,len)];
    Time_EKF2_1=[Time_EKF2_1 Runtime];
    if i==i0
        x_EKF2_1=x_RMSE(1,:);
        y_EKF2_1=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(1,magnitude,repeat,3);
    xplotEKF2_2=[xplotEKF2_2 x_RMSE(1,len)];
    xerrorplotEKF2_2=[xerrorplotEKF2_2 x_error_self_est];
    vxplotEKF2_2=[vxplotEKF2_2 x_RMSE(4,len)];
    Time_EKF2_2=[Time_EKF2_2 Runtime];
    if i==i0
        x_EKF2_2=x_RMSE(1,:);
        y_EKF2_2=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(0,magnitude,repeat,3.5);
    xplotIEKF2_1=[xplotIEKF2_1 x_RMSE(1,len)];
    xerrorplotIEKF2_1=[xerrorplotIEKF2_1 x_error_self_est];
    vxplotIEKF2_1=[vxplotIEKF2_1 x_RMSE(4,len)];
    Time_IEKF2_1=[Time_IEKF2_1 Runtime];
    if i==i0
        x_IEKF2_1=x_RMSE(1,:);
        y_IEKF2_1=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(1,magnitude,repeat,3.5);
    xplotIEKF2_2=[xplotIEKF2_2 x_RMSE(1,len)];
    xerrorplotIEKF2_2=[xerrorplotIEKF2_2 x_error_self_est];
    vxplotIEKF2_2=[vxplotIEKF2_2 x_RMSE(4,len)];
    Time_IEKF2_2=[Time_IEKF2_2 Runtime];
    if i==i0
        x_IEKF2_2=x_RMSE(1,:);
        y_IEKF2_2=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(0,magnitude,repeat,1);
    xplotUKF_1=[xplotUKF_1 x_RMSE(1,len)];
    xerrorplotUKF_1=[xerrorplotUKF_1 x_error_self_est];
    vxplotUKF_1=[vxplotUKF_1 x_RMSE(4,len)];
    Time_UKF_1=[Time_UKF_1 Runtime];
    if i==i0
        x_UKF_1=x_RMSE(1,:);
        y_UKF_1=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(1,magnitude,repeat,1);
    xplotUKF_2=[xplotUKF_2 x_RMSE(1,len)];
    xerrorplotUKF_2=[xerrorplotUKF_2 x_error_self_est];
    vxplotUKF_2=[vxplotUKF_2 x_RMSE(4,len)];
    Time_UKF_2=[Time_UKF_2 Runtime];
    if i==i0
        x_UKF_2=x_RMSE(1,:);
        y_UKF_2=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(0,magnitude,repeat,2);
    xplotCKF_1=[xplotCKF_1 x_RMSE(1,len)];
    xerrorplotCKF_1=[xerrorplotCKF_1 x_error_self_est];
    vxplotCKF_1=[vxplotCKF_1 x_RMSE(4,len)];
    Time_CKF_1=[Time_CKF_1 Runtime];
    if i==i0
        x_CKF_1=x_RMSE(1,:);
        y_CKF_1=x_RMSE(4,:);
    end

    [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(1,magnitude,repeat,2);
    xplotCKF_2=[xplotCKF_2 x_RMSE(1,len)];
    xerrorplotCKF_2=[xerrorplotCKF_2 x_error_self_est];
    vxplotCKF_2=[vxplotCKF_2 x_RMSE(4,len)];
    Time_CKF_2=[Time_CKF_2 Runtime];
    if i==i0
        x_CKF_2=x_RMSE(1,:);
        y_CKF_2=x_RMSE(4,:);
    end
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
%% plot results
f1=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
h1_2=plot(10.^scale, xplotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, xplotEKF_1, '--o', 'DisplayName', 'EKF (old)', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, xplotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, xplotEKF2_1, '--s', 'DisplayName', 'EKF2 (old)', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, xplotUKF_2, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, xplotUKF_1, '--^', 'DisplayName', 'UKF (old)', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, xplotCKF_2, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, xplotCKF_1, '--x', 'DisplayName', 'CKF (old)', 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, xplotIEKF_2, '-d', 'DisplayName', 'IEKF', 'Color', "#4DBEEE",'LineWidth',1);
h5_1=plot(10.^scale, xplotIEKF_1, '--d', 'DisplayName', 'IEKF', 'Color', "#4DBEEE",'LineWidth',1);
%h6_2=plot(10.^scale, xplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, xplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
%xlabel('Measurement noise (per unit)',FontSize=12)
ylabel('X-axis position RMSE (m)',FontSize=12)
ylim([10^-3 11])
xlim([10^-4 100])
set(gca, 'YTick', [0.001 0.01 0.1 1 10 100]);
set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
h1_1=plot(nan, nan, '--', 'Color', 'k', 'DisplayName', 'Old framework','LineWidth',1);
h2_1=plot(nan, nan, 'Color', "k", 'DisplayName', 'New framework','LineWidth',1);
leg=legend([h1_1 h2_1], 'Location','southeast',FontSize=11);
title(leg,'Line styles')
grid on
nexttile
hold on
h1_2=plot(10.^scale, vxplotEKF_2, '-o', 'DisplayName', 'EKF (new)', 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, vxplotEKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, vxplotEKF2_2, '-s', 'DisplayName', 'EKF2 (new)', 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, vxplotEKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, vxplotUKF_2, '-^', 'DisplayName', 'UKF (new)', 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, vxplotUKF_1, '--^', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, vxplotCKF_2, '-x', 'DisplayName', 'CKF (new)', 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, vxplotCKF_1, '--x', 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, vxplotIEKF_2, '-d', 'DisplayName', 'IEKF', 'Color', "#4DBEEE",'LineWidth',1);
h5_1=plot(10.^scale, vxplotIEKF_1, '--d', 'DisplayName', 'IEKF (new)', 'Color', "#4DBEEE",'LineWidth',1);
%h6_2=plot(10.^scale, vxplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, vxplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
h1_1=plot(nan, nan, 'o', 'Color', '#0072BD', 'DisplayName', 'EKF','LineWidth',1);
h2_1=plot(nan, nan, 's', 'Color', "#D95319", 'DisplayName', 'EKF2','LineWidth',1);
h3_1=plot(nan, nan, '^', 'Color', "#EDB120", 'DisplayName', 'UKF','LineWidth',1);
h4_1=plot(nan, nan, 'x', 'Color', "#7E2F8E", 'DisplayName', 'CKF','LineWidth',1);
h5_1=plot(nan, nan, 'd', 'Color', "#4DBEEE", 'DisplayName', 'IEKF','LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
xlabel('Measurement standard deviations (m)',FontSize=12)
ylabel('X-axis speed RMSE (m/s)',FontSize=12)
ylim([10^-4 1.5])
xlim([10^-4 100])
set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10 100]);
set(gca, 'YTick', [0.0001 0.001 0.01 0.1 1 10]);
leg=legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast','NumColumns',2,FontSize=11);
title(leg,'Line colors')
grid on
%exportgraphics(f1,'3D-tracking.png','Resolution',900)
%%
figure 
f4=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
h1_2=plot(0:1:len-1, x_EKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(0:1:len-1, x_EKF_1, '--o', 'DisplayName', 'EKF (old)', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(0:1:len-1, x_EKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(0:1:len-1, x_EKF2_1, '--s', 'DisplayName', 'EKF2 (old)', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(0:1:len-1, x_UKF_2, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(0:1:len-1, x_UKF_1, '--^', 'DisplayName', 'UKF (old)', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(0:1:len-1, x_CKF_2, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(0:1:len-1, x_CKF_1, '--x', 'DisplayName', 'CKF (old)', 'Color', "#7E2F8E",'LineWidth',1);
h5_1=plot(0:1:len-1, x_IEKF_1, '--d', 'DisplayName', 'IEKF', 'Color', "#4DBEEE",'LineWidth',1);
set(gca, 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
ylabel('X-axis position RMSE (m)',FontSize=12)
ylim([10^-3 25])
set(gca, 'YTick', [0.001 0.01 0.1 1 10]);
%set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
h1_1=plot(nan, nan, '--', 'Color', 'k', 'DisplayName', 'Old framework','LineWidth',1);
h2_1=plot(nan, nan, 'Color', "k", 'DisplayName', 'New framework','LineWidth',1);
leg=legend([h1_1 h2_1], 'Location','southwest',FontSize=11);
title(leg,'Line styles')
grid on
nexttile
hold on
h1_2=plot(0:1:len-1, y_EKF_2, '-o', 'DisplayName', 'EKF (new)', 'Color', '#0072BD','LineWidth',1);
h1_1=plot(0:1:len-1, y_EKF_1, '--o', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(0:1:len-1, y_EKF2_2, '-s', 'DisplayName', 'EKF2 (new)', 'Color', "#D95319",'LineWidth',1);
h2_1=plot(0:1:len-1, y_EKF2_1, '--s', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(0:1:len-1, y_UKF_2, '-^', 'DisplayName', 'UKF (new)', 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(0:1:len-1, y_UKF_1, '--^', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(0:1:len-1, y_CKF_2, '-x', 'DisplayName', 'CKF (new)', 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(0:1:len-1, y_CKF_1, '--x', 'Color', "#7E2F8E",'LineWidth',1);
h5_1=plot(0:1:len-1, y_IEKF_1, '--d', 'DisplayName', 'IEKF (new)', 'Color', "#4DBEEE",'LineWidth',1);
h1_1=plot(nan, nan, 'o', 'Color', '#0072BD', 'DisplayName', 'EKF','LineWidth',1);
h2_1=plot(nan, nan, 's', 'Color', "#D95319", 'DisplayName', 'EKF2','LineWidth',1);
h3_1=plot(nan, nan, '^', 'Color', "#EDB120", 'DisplayName', 'UKF','LineWidth',1);
h4_1=plot(nan, nan, 'x', 'Color', "#7E2F8E", 'DisplayName', 'CKF','LineWidth',1);
h5_1=plot(nan, nan, 'd', 'Color', "#4DBEEE", 'DisplayName', 'IEKF','LineWidth',1);
set(gca, 'YScale', 'log','FontSize', 12)
xlabel('Iterations',FontSize=12)
ylabel('X-axis speed RMSE (m/s)',FontSize=12)
ylim([10^-4 1.7])
%set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10 100]);
set(gca, 'YTick', [0.0001 0.001 0.01 0.1 1 10]);
leg=legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southwest','NumColumns',5,FontSize=11);
title(leg,'Line colors')
grid on
%exportgraphics(f4,'3D-tracking-convergence.png','Resolution',900)
%%
% figure 
% f3=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
% nexttile
% hold on
% h1_2=plot(10.^scale, xplotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, xplotEKF_1, '--o', 'DisplayName', ['EKF (', sprintf('%.2f', 1000*mean(Time_EKF_1)), ' ms)'], 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(10.^scale, xplotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, xplotEKF2_1, '--s', 'DisplayName', ['EKF2 (', sprintf('%.2f', 1000*mean(Time_EKF2_1)), ' ms)'], 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(10.^scale, xplotIEKF_2, '-d', 'DisplayName', 'IEKF', 'Color', "#4DBEEE",'LineWidth',1);
% h5_1=plot(10.^scale, xplotIEKF_1, '--d', 'DisplayName', ['IEKF (', sprintf('%.2f', 1000*mean(Time_IEKF_1)), ' ms)'], 'Color', "#4DBEEE",'LineWidth',1);
% h6_2=plot(10.^scale, xplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, xplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
% %xlabel('Measurement noise (per unit)',FontSize=12)
% ylabel('X-axis position RMSE (m)',FontSize=12)
% ylim([10^-3 10])
% xlim([10^-4 100])
% set(gca, 'YTick', [0.001 0.01 0.1 1 10 100]);
% set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
% %legend([h1_1 h2_1 h3_1 h4_1], 'Location','southeast',FontSize=12)
% grid on
% nexttile
% hold on
% h1_2=plot(10.^scale, vxplotEKF_2, '-o', 'DisplayName', 'EKF', 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, vxplotEKF_1, '--o', 'DisplayName', 'EKF (old)', 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(10.^scale, vxplotEKF2_2, '-s', 'DisplayName', 'EKF2', 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, vxplotEKF2_1, '--s', 'DisplayName', 'EKF2 (old)', 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(10.^scale, vxplotIEKF_2, '-d', 'DisplayName', 'IEKF', 'Color', "#4DBEEE",'LineWidth',1);
% h5_1=plot(10.^scale, vxplotIEKF_1, '--d', 'DisplayName', 'IEKF', 'Color', "#4DBEEE",'LineWidth',1);
% h6_2=plot(10.^scale, vxplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, vxplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
% xlabel('Measurement standard deviations (m)',FontSize=12)
% ylabel('X-axis speed RMSE (m/s)',FontSize=12)
% ylim([10^-4 1.5])
% xlim([10^-4 100])
% set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10 100]);
% set(gca, 'YTick', [0.0001 0.001 0.01 0.1 1 10]);
% legend([h1_2 h2_2 h5_2 h6_2], 'Location','southeast',FontSize=12)
% grid on
% exportgraphics(f3,'3D-tracking-IKF.png','Resolution',900)
%%
figure
f2=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
title("Old framework")
h1_2=plot(10.^scale, xplotEKF_1, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, xerrorplotEKF_1, '--o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, xplotEKF2_1, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, xerrorplotEKF2_1, '--s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, xplotUKF_1, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, xerrorplotUKF_1, '--^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, xplotCKF_1, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, xerrorplotCKF_1, '--x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h5_2=plot(10.^scale, xplotIEKF_1, '-d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
h5_1=plot(10.^scale, xerrorplotIEKF_1, '--d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
%h6_2=plot(10.^scale, xplotIEKF2_1, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, xerrorplotIEKF2_1, '--v', 'Color', "#A2142F",'LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
%xlabel('Measurement noise (per unit)',FontSize=12)
ylabel('X-axis position RMSE (m)',FontSize=12)
ylim([10^-3 10.5])
xlim([10^-4 100])
set(gca, 'YTick', [0.001 0.01 0.1 1 10 100]);
set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
h1_1=plot(nan, nan, '--', 'Color', 'k', 'DisplayName', 'Estimated RMSE','LineWidth',1);
h2_1=plot(nan, nan, 'Color', "k", 'DisplayName', 'Actual error','LineWidth',1);
leg=legend([h1_1 h2_1], 'Location','southeast',FontSize=11);
title(leg,'Line styles')
grid on
nexttile
hold on
title("New framework")
h1_2=plot(10.^scale, xplotEKF_2, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, xerrorplotEKF_2, '--o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, xplotEKF2_2, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, xerrorplotEKF2_2, '--s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, xplotUKF_2, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, xerrorplotUKF_2, '--^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, xplotCKF_2, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, xerrorplotCKF_2, '--x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, xplotIEKF_2, '-d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
%h5_1=plot(10.^scale, xerrorplotIEKF_2, '--d', 'DisplayName', "IEKF", 'Color', "#4DBEEE",'LineWidth',1);
%h6_2=plot(10.^scale, xplotIEKF2_2, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, xerrorplotIEKF2_2, '--v', 'Color', "#A2142F",'LineWidth',1);
h1_1=plot(nan, nan, 'o', 'Color', '#0072BD', 'DisplayName', 'EKF','LineWidth',1);
h2_1=plot(nan, nan, 's', 'Color', "#D95319", 'DisplayName', 'EKF2','LineWidth',1);
h3_1=plot(nan, nan, '^', 'Color', "#EDB120", 'DisplayName', 'UKF','LineWidth',1);
h4_1=plot(nan, nan, 'x', 'Color', "#7E2F8E", 'DisplayName', 'CKF','LineWidth',1);
h5_1=plot(nan, nan, 'd', 'Color', "#4DBEEE", 'DisplayName', 'IEKF','LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
ylabel('X-axis position RMSE (m)',FontSize=12)
xlabel('Measurement standard deviations (m)',FontSize=12)
ylim([10^-3 10.5])
xlim([10^-4 100])
set(gca, 'YTick', [0.001 0.01 0.1 1 10 100]);
set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10]);
leg=legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast','NumColumns',2,FontSize=11);
title(leg,'Line colors')
grid on
%exportgraphics(f2,'3D-tracking-error.png','Resolution',900)
%%
f5=figure();
set(f5, 'Position', [100, 100, 400, 400]);
plot(0, 0, 'p', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');hold on;
plot(p_sensor(1,:), p_sensor(2,:), 'linewidth', 2, 'Color','k','LineStyle','--'); hold on;
plot(x_true(1,:), x_true(2,:), 'linewidth', 2,'Color','k'); hold on;
legend({'Sensor 1', 'Sensor 2', 'Target'}, 'fontsize', 12);
xlabel('X-axis position (m)', 'fontsize', 12); ylabel('Y-axis position (m)', 'fontsize', 12);
xlim([-10 50])
ylim([-10 50])
axis equal;
set(gca, 'FontSize', 12)
%exportgraphics(f5,'3D-tracking-problem.png','Resolution',900)
function [time, x_true, x_RMSE, p_sensor, Runtime, x_error_self_est] = simulation(improvement,magnitude,M,KFtype)
rng(123);
N = 30; % number of time steps
dt = 1; % time between time steps
%M = 1000; % number of Monte-Carlo runs

sig_mea_true = [1.0; 1.0]*magnitude; % true value of standard deviation of measurement noise
sig_pro = [1e-3; 1e-3; 1e-3];% standard deviation of process noise
sig_mea = sig_mea_true; % user input of standard deviation of measurement noise

sig_init = [10; 10; 10; 0.1; 0.1; 0.1]; % standard deviation of initial guess

Q = [zeros(3), zeros(3); zeros(3), diag(sig_pro.^2)]; % process noise covariance matrix
R = diag(sig_mea.^2); % measurement noise covariance matrix

F = [eye(3), eye(3)*dt; zeros(3), eye(3)]; % state transition matrix
B = eye(6);
u = zeros(6,1); % noise
x_error_self_est=0;
%% true trajectory

% sensor trajectory
p_sensor = zeros(3,N+1);
for k = 1:1:N+1
    p_sensor(1,k) = 20 + 20*cos(2*pi/30 * (k-1));
    p_sensor(2,k) = 20 + 20*sin(2*pi/30 * (k-1));
    p_sensor(3,k) = 0;
end

% true target trajectory
x_true = zeros(6,N+1); 
x_true(:,1) = [10; -10; 50; 1; 2; 0]; % initial true state
for k = 2:1:N+1
    x_true(:,k) = F*x_true(:,k-1);
end

%% Kalman filter simulation
Runtime=0;
res_x_est = zeros(6,N+1,M); % Monte-Carlo estimates
res_x_err = zeros(6,N+1,M); % Monte-Carlo estimate errors
P_diag = zeros(6,N+1); % diagonal term of error covariance matrix

% filtering
for m = 1:1:M
    % initial guess
    x_est(:,1) = x_true(:,1) + normrnd(0, sig_init);
    P = [diag(sig_init(1:3).^2), zeros(3); zeros(3), diag(sig_init(4:6).^2)];
    P_diag(:,1) = diag(P);
    for k = 2:1:N+1
        u=[0; 0; 0; sig_pro].*randn(6,1);
        tic;
        if(KFtype==0)
            %%% Prediction

            % predicted state estimate
            x_est(:,k) = F*x_est(:,k-1) + B*u;
            

            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update
            pp2 = x_est(1:3,k) - p_sensor(:,k); % predicted relative position
            pp1 = x_est(1:3,k);
            
            % obtain measurement
            p = x_true(1:3,k) - p_sensor(:,k); % true relative position
            z_true = [norm(x_true(1:3,k));
                norm(p)]; % true measurement

            z = z_true + normrnd(0, sig_mea_true); % erroneous measurement

            % predicted meausrement
            
            z_p = [norm(pp1);
                norm(pp2)]; % predicted measurement

            % measurement residual
            y = z - z_p;

            % measurement matrix
            H = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];

            % Kalman gain
            K = P*H'/(R+H*P*H');

            % updated state estimate
            x_est0=x_est(:,k);
            x_est(:,k) = x_est(:,k) + K*y;

            % updated error covariance
            if(improvement==1)
                pp2 = x_est(1:3,k) - p_sensor(:,k);
                pp1 = x_est(1:3,k);
                H2 = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                    pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
                
                temp=P-K*H2*P-P*H2'*K'+K*(H2*P*H2'+R)*K';
                if(sum(diag(temp))>sum(diag(P)))
                    x_est(:,k)=x_est0;
                else
                    P=temp;
                end
            else
                P = (eye(6) - K*H)*P;
            end
        elseif(KFtype==0.5)
            %IEKF
            x_est(:,k) = F*x_est(:,k-1) + B*u;
            

            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update
            pp2 = x_est(1:3,k) - p_sensor(:,k); % predicted relative position
            pp1 = x_est(1:3,k);
            
            % obtain measurement
            p = x_true(1:3,k) - p_sensor(:,k); % true relative position
            z_true = [norm(x_true(1:3,k));
                norm(p)]; % true measurement

            z = z_true + normrnd(0, sig_mea_true); % erroneous measurement

            % predicted meausrement
            
            z_p = [norm(pp1);
                norm(pp2)]; % predicted measurement

            % measurement residual
            y = z - z_p;

            % measurement matrix
            H = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];

            % Kalman gain
            K = P*H'/(R+H*P*H');
            % updated state estimate
            x_est0=x_est(:,k);
            x_est(:,k) = x_est(:,k) + K*y;
            %iteration
            steplennow = norm(K*y);
            iter=1;
            change=1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                pp2 = x_est(1:3,k) - p_sensor(:,k);
                pp1 = x_est(1:3,k);
                z_p = [norm(pp1);
                    norm(pp2)]; % predicted measurement
                y = z - z_p;
                H2 = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                    pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
                S = H2*P*H2.'+ R;
                K2 = P*H2.'*S^(-1);
                dx = x_est0 + K2*(y-H2*(x_est0-x_est(:,k))) - x_est(:,k);
                steplen_previous=steplennow;
                steplennow=norm(dx);
                if(steplen_previous<steplennow)
                    break;
                end
                change = max(abs(dx./x_est(:,k)));
                x_est(:,k) = x_est(:,k) + dx;
                K = K2;
                H = H2;
            end
            % updated error covariance
            if(improvement==1)
                pp2 = x_est(1:3,k) - p_sensor(:,k);
                pp1 = x_est(1:3,k);
                H2 = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                    pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
                
                temp=P-K*H2*P-P*H2'*K'+K*(H2*P*H2'+R)*K';
                if(sum(diag(temp))>sum(diag(P)))
                    x_est(:,k)=x_est0;
                else
                    P=temp;
                end
            else
                P = (eye(6) - K*H)*P;
            end
        elseif(KFtype==1)
            L=real(sqrtm(P));
            state=x_est(:,k-1);
            n=6;
            lambda=(1e-6-1)*n;
            alpha=1e-3;
            beta=2;
            mnum=2;
            states=zeros(n,n*2+1);
            states(:,1)=state;
            for i=1:n
                states(:,1+i)=state+sqrt(lambda+n)*L(:,i);
                states(:,n+i+1)=state-sqrt(lambda+n)*L(:,i);  
            end
            for i=1:2*n+1
                states(:,i)=F*states(:,i) + B*u;
                if(i==1)
                    state=states(:,i)*lambda/(n+lambda);
                else
                    state=state+states(:,i)/(n+lambda)/2;
                end
            end
            P=(states(:,1)-state)*(states(:,1)-state).'*(lambda/(n+lambda)+1-alpha^2+beta)+Q;
            for i=2:2*n+1
                P=P+(state-states(:,i))*(state-states(:,i)).'/(n+lambda)/2;
            end
            % obtain measurement
            p = x_true(1:3,k) - p_sensor(:,k); % true relative position
            z_true = [norm(x_true(1:3,k));
                norm(p)]; % true measurement

            z = z_true + normrnd(0, sig_mea_true); % erroneous measurement
            % Predict Measurement From Propagated Sigma Points
            measures=zeros(mnum,2*n+1);
            for i=1:2*n+1
                pp2 = states(1:3,i) - p_sensor(:,k); % predicted relative position
                pp1 = states(1:3,i);
                measures(:,i) = [norm(pp1);
                    norm(pp2)]; % predicted measurement
                if(i==1)
                    m_exp=lambda/(n+lambda)*measures(:,i);
                else
                    m_exp=m_exp+1/(n+lambda)/2*measures(:,i);
                end
            end
            % Estimate Mean And Covariance of Predicted Measurement
            Py=(lambda/(n+lambda)+1-alpha^2+beta)*(measures(:,1)-m_exp)*(measures(:,1)-m_exp).'+R;
            Pxy=(lambda/(n+lambda)+1-alpha^2+beta)*(states(:,1)-state)*(measures(:,1)-m_exp).';
            for i=2:2*n+1
                Py=Py+1/(n+lambda)/2*(measures(:,i)-m_exp)*(measures(:,i)-m_exp).';
                Pxy=Pxy+1/(n+lambda)/2*(states(:,i)-state)*(measures(:,i)-m_exp).';
            end
            %kalman gain
            K=Pxy/Py;
            %update
            state0=state;
            state=state+K*(z-m_exp);
            if(improvement==1)
                L=real(sqrtm(P));
                states=zeros(n,n*2+1);
                states(:,1)=state;
                for i=1:n
                    states(:,1+i)=state+sqrt(lambda+n)*L(:,i);
                    states(:,n+i+1)=state-sqrt(lambda+n)*L(:,i);
                end
                measures=zeros(mnum,2*n+1);
                for i=1:2*n+1
                    pp2 = states(1:3,i) - p_sensor(:,k); % predicted relative position
                    pp1 = states(1:3,i);
                    measures(:,i) = [norm(pp1);
                        norm(pp2)]; % predicted measurement
                    if(i==1)
                        m_exp=lambda/(n+lambda)*measures(:,i);
                    else
                        m_exp=m_exp+1/(n+lambda)/2*measures(:,i);
                    end
                end
                % Estimate Mean And Covariance of Predicted Measurement
                Py=(lambda/(n+lambda)+1-alpha^2+beta)*(measures(:,1)-m_exp)*(measures(:,1)-m_exp).'+R;
                Pxy=(lambda/(n+lambda)+1-alpha^2+beta)*(states(:,1)-state)*(measures(:,1)-m_exp).';
                for i=2:2*n+1
                    Py=Py+1/(n+lambda)/2*(measures(:,i)-m_exp)*(measures(:,i)-m_exp).';
                    Pxy=Pxy+1/(n+lambda)/2*(states(:,i)-state)*(measures(:,i)-m_exp).';
                end
                temp=P+K*Py*K.'-Pxy*K.'-K*Pxy.';
                if(sum(diag(temp))<sum(diag(P)))
                    P=temp;
                else
                    state=state0;
                end
            else
                P=P-K*Py*K.';
            end
            x_est(:,k)=state;
        elseif(KFtype==2)
            %CKF
            L=real(sqrtm(P));
            state=x_est(:,k-1);
            n=6;
            mnum=2;
            states=zeros(n,n*2);
            for i=1:n
                states(:,i)=state+sqrt(n)*L(:,i);
                states(:,n+i)=state-sqrt(n)*L(:,i);  
            end
            state=0;
            for i=1:2*n
                states(:,i)=F*states(:,i) + B*u;
                state=state+states(:,i)/n/2;
            end
            P=Q;
            for i=1:2*n
                P=P+(state-states(:,i))*(state-states(:,i)).'/n/2;
            end
            % obtain measurement
            p = x_true(1:3,k) - p_sensor(:,k); % true relative position
            z_true = [norm(x_true(1:3,k));
                norm(p)]; % true measurement

            z = z_true + normrnd(0, sig_mea_true); % erroneous measurement
            %update sigma points
            L=real(sqrtm(P));
            states=zeros(n,n*2);
            for i=1:n
                states(:,i)=state+sqrt(n)*L(:,i);
                states(:,n+i)=state-sqrt(n)*L(:,i);  
            end
            % Predict Measurement From Propagated Sigma Points
            measures=zeros(mnum,2*n);
            m_exp=0;
            for i=1:2*n
                pp2 = states(1:3,i) - p_sensor(:,k); % predicted relative position
                pp1 = states(1:3,i);
                measures(:,i) = [norm(pp1);
                    norm(pp2)]; % predicted measurement
                m_exp=m_exp+1/n/2*measures(:,i);
            end
            % Estimate Mean And Covariance of Predicted Measurement
            Py=R;
            Pxy=0;
            for i=1:2*n
                Py=Py+1/n/2*(measures(:,i)-m_exp)*(measures(:,i)-m_exp).';
                Pxy=Pxy+1/n/2*(states(:,i)-state)*(measures(:,i)-m_exp).';
            end
            %kalman gain
            K=Pxy/Py;
            %update
            state0=state;
            state=state+K*(z-m_exp);
            if(improvement==1)
                states=zeros(n,n*2);
                for i=1:n
                    states(:,i)=state+sqrt(n)*L(:,i);
                    states(:,n+i)=state-sqrt(n)*L(:,i);
                end
                % Predict Measurement From Propagated Sigma Points
                measures=zeros(mnum,2*n);
                m_exp=0;
                for i=1:2*n
                    pp2 = states(1:3,i) - p_sensor(:,k); % predicted relative position
                    pp1 = states(1:3,i);
                    measures(:,i) = [norm(pp1);
                        norm(pp2)]; % predicted measurement
                    m_exp=m_exp+1/n/2*measures(:,i);
                end
                % Estimate Mean And Covariance of Predicted Measurement
                Py=R;
                Pxy=0;
                for i=1:2*n
                    Py=Py+1/n/2*(measures(:,i)-m_exp)*(measures(:,i)-m_exp).';
                    Pxy=Pxy+1/n/2*(states(:,i)-state)*(measures(:,i)-m_exp).';
                end
                
                temp=P+K*Py*K.'-Pxy*K.'-K*Pxy.';
                if(sum(diag(temp))<sum(diag(P)))
                    P=temp;
                else
                    state=state0;
                end
            else
                P=P-K*Py*K.';
            end
            x_est(:,k)=state;
        elseif(KFtype==3)
            %2-EKF
            %%% Prediction

            % predicted state estimate
            x_est(:,k) = F*x_est(:,k-1) + B*u;
            n=length(x_est(:,k));

            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update
            pp2 = x_est(1:3,k) - p_sensor(:,k); % predicted relative position
            pp1 = x_est(1:3,k);
            
            % obtain measurement
            
            p = x_true(1:3,k) - p_sensor(:,k); % true relative position
            z_true = [norm(x_true(1:3,k));
                norm(p)]; % true measurement

            z = z_true + normrnd(0, sig_mea_true); % erroneous measurement

            % predicted meausrement
            mnum=2;
            z_p = [norm(pp1);
                norm(pp2)]; % predicted measurement

            % measurement matrix
            H = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
            Hxx=zeros(n,n,mnum);
            Hxx(1:3,1:3,1)=[1/norm(pp1)-pp1(1)^2/norm(pp1)^3 -pp1(1)*pp1(2)/norm(pp1)^3 -pp1(1)*pp1(3)/norm(pp1)^3;
               -pp1(1)*pp1(2)/norm(pp1)^3 1/norm(pp1)-pp1(2)^2/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3;
               -pp1(1)*pp1(3)/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3 1/norm(pp1)-pp1(3)^2/norm(pp1)^3];
            Hxx(1:3,1:3,2)=[1/norm(pp2)-pp2(1)^2/norm(pp2)^3 -pp2(1)*pp2(2)/norm(pp2)^3 -pp2(1)*pp2(3)/norm(pp2)^3;
               -pp2(1)*pp2(2)/norm(pp2)^3 1/norm(pp2)-pp2(2)^2/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3;
               -pp2(1)*pp2(3)/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3 1/norm(pp2)-pp2(3)^2/norm(pp2)^3];
            S=H*P*H.'+R;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*P*Hxx(:,:,jj)*P));
                end
            end
            % Kalman gain
            K = P*H.'*S^(-1);
            % measurement residual
            residual=z-z_p;
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*P));
            end
            % updated state estimate
            x_est0 = x_est(:,k);
            x_est(:,k) = x_est(:,k) + K*residual;

            % updated error covariance
            if(improvement==1)
                pp2 = x_est(1:3,k) - p_sensor(:,k);
                pp1 = x_est(1:3,k);
                H2 = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                    pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
                Hxx2=zeros(n,n,mnum);
                Hxx2(1:3,1:3,1)=[1/norm(pp1)-pp1(1)^2/norm(pp1)^3 -pp1(1)*pp1(2)/norm(pp1)^3 -pp1(1)*pp1(3)/norm(pp1)^3;
               -pp1(1)*pp1(2)/norm(pp1)^3 1/norm(pp1)-pp1(2)^2/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3;
               -pp1(1)*pp1(3)/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3 1/norm(pp1)-pp1(3)^2/norm(pp1)^3];
                Hxx2(1:3,1:3,2)=[1/norm(pp2)-pp2(1)^2/norm(pp2)^3 -pp2(1)*pp2(2)/norm(pp2)^3 -pp2(1)*pp2(3)/norm(pp2)^3;
                   -pp2(1)*pp2(2)/norm(pp2)^3 1/norm(pp2)-pp2(2)^2/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3;
                   -pp2(1)*pp2(3)/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3 1/norm(pp2)-pp2(3)^2/norm(pp2)^3];
                S2=H2*P*H2.'+R;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*P*Hxx2(:,:,jj)*P));
                    end
                end
                temp=P+K*S2*K.'-P*H2.'*K.'-K*H2*P;
                if(sum(diag(temp))>sum(diag(P)))
                    x_est(:,k)=x_est0;
                else
                    P=temp;
                end

            else
                P = P-K*S*K.';
            end
        else
            %2-EKF
            %%% Prediction

            % predicted state estimate
            x_est(:,k) = F*x_est(:,k-1) + B*u;
            n=length(x_est(:,k));

            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update
            pp2 = x_est(1:3,k) - p_sensor(:,k); % predicted relative position
            pp1 = x_est(1:3,k);

            % obtain measurement

            p = x_true(1:3,k) - p_sensor(:,k); % true relative position
            z_true = [norm(x_true(1:3,k));
                norm(p)]; % true measurement

            z = z_true + normrnd(0, sig_mea_true); % erroneous measurement

            % predicted meausrement
            mnum=2;
            z_p = [norm(pp1);
                norm(pp2)]; % predicted measurement

            % measurement matrix
            H = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
            Hxx=zeros(n,n,mnum);
            Hxx(1:3,1:3,1)=[1/norm(pp1)-pp1(1)^2/norm(pp1)^3 -pp1(1)*pp1(2)/norm(pp1)^3 -pp1(1)*pp1(3)/norm(pp1)^3;
                -pp1(1)*pp1(2)/norm(pp1)^3 1/norm(pp1)-pp1(2)^2/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3;
                -pp1(1)*pp1(3)/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3 1/norm(pp1)-pp1(3)^2/norm(pp1)^3];
            Hxx(1:3,1:3,2)=[1/norm(pp2)-pp2(1)^2/norm(pp2)^3 -pp2(1)*pp2(2)/norm(pp2)^3 -pp2(1)*pp2(3)/norm(pp2)^3;
                -pp2(1)*pp2(2)/norm(pp2)^3 1/norm(pp2)-pp2(2)^2/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3;
                -pp2(1)*pp2(3)/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3 1/norm(pp2)-pp2(3)^2/norm(pp2)^3];
            S=H*P*H.'+R;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*P*Hxx(:,:,jj)*P));
                end
            end
            % Kalman gain
            K = P*H.'*S^(-1);
            % measurement residual
            residual=z-z_p;
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*P));
            end
            % updated state estimate
            x_est0 = x_est(:,k);
            x_est(:,k) = x_est(:,k) + K*residual;
            
            % iteration
            steplennow = norm(K*residual);
            iter=1;
            change=1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                pp2 = x_est(1:3,k) - p_sensor(:,k);
                pp1 = x_est(1:3,k);
                z_p = [norm(pp1);
                    norm(pp2)]; % predicted measurement
                H2 = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                    pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
                Hxx2=zeros(n,n,mnum);
                Hxx2(1:3,1:3,1)=[1/norm(pp1)-pp1(1)^2/norm(pp1)^3 -pp1(1)*pp1(2)/norm(pp1)^3 -pp1(1)*pp1(3)/norm(pp1)^3;
                    -pp1(1)*pp1(2)/norm(pp1)^3 1/norm(pp1)-pp1(2)^2/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3;
                    -pp1(1)*pp1(3)/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3 1/norm(pp1)-pp1(3)^2/norm(pp1)^3];
                Hxx2(1:3,1:3,2)=[1/norm(pp2)-pp2(1)^2/norm(pp2)^3 -pp2(1)*pp2(2)/norm(pp2)^3 -pp2(1)*pp2(3)/norm(pp2)^3;
                    -pp2(1)*pp2(2)/norm(pp2)^3 1/norm(pp2)-pp2(2)^2/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3;
                    -pp2(1)*pp2(3)/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3 1/norm(pp2)-pp2(3)^2/norm(pp2)^3];
                residual=z-z_p;
                for ii=1:mnum
                    residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*P));
                    residual(ii)=residual(ii)-0.5*(x_est0-x_est(:,k)).'*Hxx2(:,:,ii)*(x_est0-x_est(:,k));
                end
                S2 = H2*P*H2.'+ R;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*P*Hxx2(:,:,jj)*P));
                    end
                end
                K2 = P*H2.'*S2^(-1);
                dx = x_est0 + K2*(residual-H2*(x_est0-x_est(:,k))) - x_est(:,k);
                steplen_previous=steplennow;
                steplennow=norm(dx);
                if(steplen_previous<steplennow)
                    break;
                end
                change = max(abs(dx./x_est(:,k)));
                x_est(:,k) = x_est(:,k) + dx;
                K = K2;
                S = S2;
            end
            % updated error covariance
            if(improvement==1)
                pp2 = x_est(1:3,k) - p_sensor(:,k);
                pp1 = x_est(1:3,k);
                H2 = [pp1(1)/norm(pp1), pp1(2)/norm(pp1), pp1(3)/norm(pp1), zeros(1,3);
                    pp2(1)/norm(pp2), pp2(2)/norm(pp2), pp2(3)/norm(pp2), zeros(1,3)];
                Hxx2=zeros(n,n,mnum);
                Hxx2(1:3,1:3,1)=[1/norm(pp1)-pp1(1)^2/norm(pp1)^3 -pp1(1)*pp1(2)/norm(pp1)^3 -pp1(1)*pp1(3)/norm(pp1)^3;
                    -pp1(1)*pp1(2)/norm(pp1)^3 1/norm(pp1)-pp1(2)^2/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3;
                    -pp1(1)*pp1(3)/norm(pp1)^3 -pp1(2)*pp1(3)/norm(pp1)^3 1/norm(pp1)-pp1(3)^2/norm(pp1)^3];
                Hxx2(1:3,1:3,2)=[1/norm(pp2)-pp2(1)^2/norm(pp2)^3 -pp2(1)*pp2(2)/norm(pp2)^3 -pp2(1)*pp2(3)/norm(pp2)^3;
                    -pp2(1)*pp2(2)/norm(pp2)^3 1/norm(pp2)-pp2(2)^2/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3;
                    -pp2(1)*pp2(3)/norm(pp2)^3 -pp2(2)*pp2(3)/norm(pp2)^3 1/norm(pp2)-pp2(3)^2/norm(pp2)^3];
                S2=H2*P*H2.'+R;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*P*Hxx2(:,:,jj)*P));
                    end
                end
                temp=P+K*S2*K.'-P*H2.'*K.'-K*H2*P;
                if(sum(diag(temp))>sum(diag(P)))
                    x_est(:,k)=x_est0;
                else
                    P=temp;
                end

            else
                P = P-K*S*K.';
            end

        end
        temp=toc;
        Runtime = Runtime + temp/M;
        P_diag(:,k) = diag(P);
    end
    
    res_x_est(:,:,m) = x_est;
    res_x_err(:,:,m) = x_est - x_true;
    x_error_self_est = x_error_self_est + P(1,1)/M;
    time = (0:1:N)*dt;
end
x_error_self_est = sqrt(x_error_self_est);
%% get result statistics

x_RMSE = zeros(6,N+1); % root mean square error
for k = 1:1:N+1
    x_RMSE(1,k) = sqrt(mean(res_x_err(1,k,:).^2,3));
    x_RMSE(2,k) = sqrt(mean(res_x_err(2,k,:).^2,3));
    x_RMSE(3,k) = sqrt(mean(res_x_err(3,k,:).^2,3));
    x_RMSE(4,k) = sqrt(mean(res_x_err(4,k,:).^2,3));
    x_RMSE(5,k) = sqrt(mean(res_x_err(5,k,:).^2,3));
    x_RMSE(6,k) = sqrt(mean(res_x_err(6,k,:).^2,3));
end
end



