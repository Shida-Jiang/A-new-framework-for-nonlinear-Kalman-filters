%%% This code is adapted from:
% Y. Kim and H. Bang, Introduction to Kalman Filter and Its Applications, InTechOpen, 2018

% Application 2:  Terrain-referenced navigation

close all
clc
clear

repeat=10000;%The number of times that the simulation is repeated at each measurement noise setup
%This simulation takes about 30 minutes, you can reduce the repeat times to reduce the run time
%% settings
scale=-1:0.5:4; % The range of measurement noise is 10^scale
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
    magnitude=10^(i-3);
    [~,x_RMSE,Runtime]=simulation(0,magnitude,repeat,0);
    len=length(x_RMSE(1,:));
    xplotEKF_1=[xplotEKF_1 x_RMSE(1,len)];
    yplotEKF_1=[yplotEKF_1 x_RMSE(2,len)];
    Time_EKF_1=[Time_EKF_1 Runtime];

    [~,x_RMSE,Runtime]=simulation(1,magnitude,repeat,0);
    len=length(x_RMSE(1,:));
    xplotEKF_2=[xplotEKF_2 x_RMSE(1,len)];
    yplotEKF_2=[yplotEKF_2 x_RMSE(2,len)];
    Time_EKF_2=[Time_EKF_2 Runtime];

    [~,x_RMSE,Runtime]=simulation(0,magnitude,repeat,0.5);
    len=length(x_RMSE(1,:));
    xplotIEKF_1=[xplotIEKF_1 x_RMSE(1,len)];
    yplotIEKF_1=[yplotIEKF_1 x_RMSE(2,len)];
    Time_IEKF_1=[Time_IEKF_1 Runtime];

    [~,x_RMSE,Runtime]=simulation(1,magnitude,repeat,0.5);
    len=length(x_RMSE(1,:));
    xplotIEKF_2=[xplotIEKF_2 x_RMSE(1,len)];
    yplotIEKF_2=[yplotIEKF_2 x_RMSE(2,len)];
    Time_IEKF_2=[Time_IEKF_2 Runtime];

    [~,x_RMSE,Runtime]=simulation(0,magnitude,repeat,3);
    xplotEKF2_1=[xplotEKF2_1 x_RMSE(1,len)];
    yplotEKF2_1=[yplotEKF2_1 x_RMSE(2,len)];
    Time_EKF2_1=[Time_EKF2_1 Runtime];

    [~,x_RMSE,Runtime]=simulation(1,magnitude,repeat,3);
    xplotEKF2_2=[xplotEKF2_2 x_RMSE(1,len)];
    yplotEKF2_2=[yplotEKF2_2 x_RMSE(2,len)];
    Time_EKF2_2=[Time_EKF2_2 Runtime];

    [~,x_RMSE,Runtime]=simulation(0,magnitude,repeat,3.5);
    xplotIEKF2_1=[xplotIEKF2_1 x_RMSE(1,len)];
    yplotIEKF2_1=[yplotIEKF2_1 x_RMSE(2,len)];
    Time_IEKF2_1=[Time_IEKF2_1 Runtime];

    [~,x_RMSE,Runtime]=simulation(1,magnitude,repeat,3.5);
    xplotIEKF2_2=[xplotIEKF2_2 x_RMSE(1,len)];
    yplotIEKF2_2=[yplotIEKF2_2 x_RMSE(2,len)];
    Time_IEKF2_2=[Time_IEKF2_2 Runtime];

    [~,x_RMSE,Runtime]=simulation(0,magnitude,repeat,1);
    xplotUKF_1=[xplotUKF_1 x_RMSE(1,len)];
    yplotUKF_1=[yplotUKF_1 x_RMSE(2,len)];
    Time_UKF_1=[Time_UKF_1 Runtime];

    [~,x_RMSE,Runtime]=simulation(1,magnitude,repeat,1);
    xplotUKF_2=[xplotUKF_2 x_RMSE(1,len)];
    yplotUKF_2=[yplotUKF_2 x_RMSE(2,len)];
    Time_UKF_2=[Time_UKF_2 Runtime];

    [~,x_RMSE,Runtime]=simulation(0,magnitude,repeat,2);
    xplotCKF_1=[xplotCKF_1 x_RMSE(1,len)];
    yplotCKF_1=[yplotCKF_1 x_RMSE(2,len)];
    Time_CKF_1=[Time_CKF_1 Runtime];

    [time,x_RMSE,Runtime]=simulation(1,magnitude,repeat,2);
    xplotCKF_2=[xplotCKF_2 x_RMSE(1,len)];
    yplotCKF_2=[yplotCKF_2 x_RMSE(2,len)];
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

%% plot
x = 0:10:1200;
y = 0:10:1200;

% Create a grid of x and y values
[X, Y] = meshgrid(x, y);

% Define the plane equation
Z = 100*sin(sqrt((X/40).^2 + (Y/40).^2));

% Plot the plane using surf function
f2=figure;
surf(X, Y, Z);

% Add colorbar to represent z values
colorbar;

% Add labels and title
xlabel('X-axis (km)');
ylabel('Y-axis (km)');
zlabel('Z-axis');
title("Altitude map")
view(2)
%exportgraphics(f2,'Terrain-referenced-navigation-x.png','Resolution',900)
%%
f1=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
h1_2=plot(10.^scale, xplotEKF_2*1000, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, xplotEKF_1*1000, '--o', 'DisplayName', "EKF (old)", 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, xplotEKF2_2*1000, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, xplotEKF2_1*1000, '--s', 'DisplayName', "EKF2 (old)", 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, xplotUKF_2*1000, '-^', 'DisplayName', "UKF", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, xplotUKF_1*1000, '--^', 'DisplayName', "UKF (old)", 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, xplotCKF_2*1000, '-x', 'DisplayName', "CKF", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, xplotCKF_1*1000, '--x', 'DisplayName', "CKF (old)", 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, xplotIEKF_2*1000, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
h5_1=plot(10.^scale, xplotIEKF_1*1000, '--d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
%h6_2=plot(10.^scale, xplotIEKF2_2*1000, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, xplotIEKF2_1*1000, '--v', 'Color', "#A2142F",'LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
%xlabel('Measurement noise (per unit)',FontSize=11)
ylabel('Abscissa RMSE (m)',FontSize=13)
%legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast',FontSize=11)
ylim([0.01 1000])
set(gca, 'YTick', [0.01 0.1 1 10 100 1000]);
xlim([10^-1 10^4])
grid on
nexttile
hold on
h1_2=plot(10.^scale, yplotEKF_2*1000, '-o', 'DisplayName', "EKF (new)", 'Color', '#0072BD','LineWidth',1);
h1_1=plot(10.^scale, yplotEKF_1*1000, '--o', 'Color', '#0072BD','LineWidth',1);
h2_2=plot(10.^scale, yplotEKF2_2*1000, '-s', 'DisplayName', "EKF2 (new)", 'Color', "#D95319",'LineWidth',1);
h2_1=plot(10.^scale, yplotEKF2_1*1000, '--s', 'Color', "#D95319",'LineWidth',1);
h3_2=plot(10.^scale, yplotUKF_2*1000, '-^', 'DisplayName', "UKF (new)", 'Color', "#EDB120",'LineWidth',1);
h3_1=plot(10.^scale, yplotUKF_1*1000, '--^', 'Color', "#EDB120",'LineWidth',1);
h4_2=plot(10.^scale, yplotCKF_2*1000, '-x', 'DisplayName', "CKF (new)", 'Color', "#7E2F8E",'LineWidth',1);
h4_1=plot(10.^scale, yplotCKF_1*1000, '--x', 'Color', "#7E2F8E",'LineWidth',1);
%h5_2=plot(10.^scale, yplotIEKF_2*1000, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
h5_1=plot(10.^scale, yplotIEKF_1*1000, '--d', 'Color', '#4DBEEE','LineWidth',1);
%h6_2=plot(10.^scale, yplotIEKF2_2*1000, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
%h6_1=plot(10.^scale, yplotIEKF2_1*1000, '--v', 'Color', "#A2142F",'LineWidth',1);
h1_1=plot(nan, nan, 'o', 'Color', '#0072BD', 'DisplayName', 'EKF','LineWidth',1);
h2_1=plot(nan, nan, 's', 'Color', "#D95319", 'DisplayName', 'EKF2','LineWidth',1);
h3_1=plot(nan, nan, '^', 'Color', "#EDB120", 'DisplayName', 'UKF','LineWidth',1);
h4_1=plot(nan, nan, 'x', 'Color', "#7E2F8E", 'DisplayName', 'CKF','LineWidth',1);
h5_1=plot(nan, nan, 'd', 'Color', "#4DBEEE", 'DisplayName', 'IEKF','LineWidth',1);
set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
xlabel('Measurement standard deviation (m)',FontSize=13)
ylabel('Ordinate RMSE (m)',FontSize=13)
legend([h1_1 h2_1 h3_1 h4_1 h5_1], 'Location','southeast',FontSize=11)
ylim([0.1 2000])
xlim([10^-1 10^4])
set(gca, 'YTick', [0.1 1 10 100 1000]);
grid on
% exportgraphics(f1,'Terrain-referenced-navigation.png','Resolution',900)

%%
% figure
% f2=tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
% nexttile
% hold on
% h1_2=plot(10.^scale, xplotEKF_2*1000, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, xplotEKF_1*1000, '--o', 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(10.^scale, xplotEKF2_2*1000, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, xplotEKF2_1*1000, '--s', 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(10.^scale, xplotIEKF_2*1000, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
% h5_1=plot(10.^scale, xplotIEKF_1*1000, '--d', 'Color', '#4DBEEE','LineWidth',1);
% h6_2=plot(10.^scale, xplotIEKF2_2*1000, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, xplotIEKF2_1*1000, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log', 'XTickLabel', [],'FontSize', 12)
% %xlabel('Measurement noise (per unit)',FontSize=11)
% ylabel('Abscissa RMSE (m)',FontSize=13)
% %legend([h1_2 h2_2 h3_2 h4_2], 'Location','southeast',FontSize=12)
% ylim([1 2000])
% set(gca, 'YTick', [1 10 100 1000]);
% xlim([10^-1 10^4])
% grid on
% nexttile
% hold on
% h1_2=plot(10.^scale, yplotEKF_2*1000, '-o', 'DisplayName', "EKF", 'Color', '#0072BD','LineWidth',1);
% h1_1=plot(10.^scale, yplotEKF_1*1000, '--o', 'Color', '#0072BD','LineWidth',1);
% h2_2=plot(10.^scale, yplotEKF2_2*1000, '-s', 'DisplayName', "EKF2", 'Color', "#D95319",'LineWidth',1);
% h2_1=plot(10.^scale, yplotEKF2_1*1000, '--s', 'Color', "#D95319",'LineWidth',1);
% h5_2=plot(10.^scale, yplotIEKF_2*1000, '-d', 'DisplayName', "IEKF", 'Color', '#4DBEEE','LineWidth',1);
% h5_1=plot(10.^scale, yplotIEKF_1*1000, '--d', 'Color', '#4DBEEE','LineWidth',1);
% h6_2=plot(10.^scale, yplotIEKF2_2*1000, '-v', 'DisplayName', "IEKF2", 'Color', "#A2142F",'LineWidth',1);
% h6_1=plot(10.^scale, yplotIEKF2_1*1000, '--v', 'Color', "#A2142F",'LineWidth',1);
% set(gca, 'XScale', 'log', 'YScale', 'log','FontSize', 12)
% xlabel('Measurement standard deviation (m)',FontSize=13)
% ylabel('Ordinate RMSE (m)',FontSize=13)
% legend([h1_2 h2_2 h5_2 h6_2], 'Location','southeast',FontSize=12)
% ylim([0.1 2000])
% xlim([10^-1 10^4])
% set(gca, 'YTick', [0.1 1 10 100 1000]);
% grid on
% exportgraphics(f2,'Terrain-referenced-navigation-IKF.png','Resolution',900)


%% function to get DEM data

function h = height(x,y)
h=1000*sin(sqrt((x/40).^2 + (y/40).^2));
end

function H = dheight(x,y)
H=1000*cos(sqrt((x/40)^2 + (y/40)^2))/sqrt((x/40)^2 + (y/40)^2)*[x/1600 y/1600];
end

function H = ddheight(x,y)
t=sqrt((x/40)^2 + (y/40)^2);
H=10*[cos(t)/16/t+x/16*(-sin(t)/t-cos(t)/t^2)*x/1600/t x/16*(-sin(t)/t-cos(t)/t^2)*y/1600/t;
    x/16*(-sin(t)/t-cos(t)/t^2)*y/1600/t cos(t)/16/t+y/16*(-sin(t)/t-cos(t)/t^2)*y/1600/t];
end
%% main simulation function
function [time,x_RMSE,Runtime]=simulation(improvement,magnitude,M,KFtype)
rng(123)
N = 100; % number of time steps
dt = 1; % time between time steps
% M is the number of Monte-Carlo runs
sig_pro_true = 0.5*[1e-3; 1e-3]; % true value of standard deviation of process noise
sig_mea_true = magnitude; % true value of standard deviation of measurement noise

sig_pro = sig_pro_true; % user input of standard deviation of process noise
sig_mea = sig_mea_true; % user input of standard deviation of measurement noise

sig_init = [1; 1]; % standard deviation of initial guess

Q = diag(sig_pro.^2); % process noise covariance matrix
R = diag(sig_mea.^2); % measurement noise covariance matrix

F = eye(2); % state transition matrix
B = eye(2); % control-input matrix

%% true trajectory

% aircraft trajectory
x_true = zeros(2,N+1); 
x_true(:,1) = [10; 10]; % initial true state
u = [0.5; 0];
for k = 2:1:N+1
    x_true(:,k) = F*x_true(:,k-1) + B*u;
end

%% Kalman filter simulation
Runtime=0;
res_x_est = zeros(2,N+1,M); % Monte-Carlo estimates
res_x_err = zeros(2,N+1,M); % Monte-Carlo estimate errors
P_diag = zeros(2,N+1); % diagonal term of error covariance matrix

% filtering
for m = 1:1:M
    % initial guess
    x_est(:,1) = x_true(:,1) + normrnd(0, sig_init);
    P = diag(sig_init.^2);
    P_diag(:,1) = diag(P);
    for k = 2:1:N+1
        
        %%% Prediction
        u_p = u + normrnd(0, sig_pro_true);
        % obtain measurement
        z = height(x_true(1,k), x_true(2,k)) + normrnd(0, sig_mea_true);
        tic;
        if(KFtype==0)
            % translation
            

            % predicted state estimate
            x_est(:,k) = F*x_est(:,k-1) + B*u_p;

            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update

            % predicted meausrement
            z_p = height(x_est(1,k), x_est(2,k));

            % measurement residual
            y = z - z_p;

            % measurement matrix
            H = dheight(x_est(1,k),x_est(2,k));

            % Kalman gain
            K = P*H'/(R+H*P*H');

            % updated state estimate
            x_est0 = x_est(:,k);
            x_est(:,k) = x_est(:,k) + K*y;

            % updated error covariance
            if(improvement==1)
                H2 = dheight(x_est(1,k),x_est(2,k));
                temp=P-K*H2*P-P*H2'*K'+K*(H2*P*H2'+R)*K';
                if(sum(diag(temp))>sum(diag(P)))
                    x_est(:,k)=x_est0;
                else
                    P=temp;
                end
            else
                P = (eye(2) - K*H)*P;
            end
        elseif(KFtype==0.5)
            % translation
            

            % predicted state estimate
            x_est(:,k) = F*x_est(:,k-1) + B*u_p;

            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update

            % predicted meausrement
            z_p = height(x_est(1,k), x_est(2,k));

            % measurement residual
            y = z - z_p;

            % measurement matrix
            H = dheight(x_est(1,k),x_est(2,k));

            % Kalman gain
            K = P*H'/(R+H*P*H');

            % updated state estimate
            x_est0 = x_est(:,k);
            x_est(:,k) = x_est(:,k) + K*y;
            %IEKF
            steplennow = norm(K*y);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                z_p = height(x_est(1,k), x_est(2,k));
                y = z - z_p;
                H2 = dheight(x_est(1,k),x_est(2,k));
                S = H2*P*H2.'+ R;
                K2 = P*H2.'*S^(-1);
                dx = x_est0 + K2*(y-H2*(x_est0-x_est(:,k))) - x_est(:,k);
                steplen_previous=steplennow;
                steplennow=norm(dx);
                if(steplen_previous<steplennow)
                    break;
                else
                    change = max(abs(dx./x_est(:,k)));
                    x_est(:,k) = x_est(:,k) + dx;
                    K = K2;
                    H = H2;
                end
            end
            % updated error covariance
            if(improvement==1)
                H2 = dheight(x_est(1,k),x_est(2,k));
                temp=P-K*H2*P-P*H2'*K'+K*(H2*P*H2'+R)*K';
                if(sum(diag(temp))>sum(diag(P)))
                    x_est(:,k)=x_est0;
                else
                    P=temp;
                end
            else
                P = (eye(2) - K*H)*P;
            end
        elseif(KFtype==1)
            L=real(sqrtm(P));
            state=x_est(:,k-1);
            n=2;
            lambda=(1e-6-1)*n;
            alpha=1e-3;
            beta=2;
            mnum=1;
            states=zeros(n,n*2+1);
            states(:,1)=state;
            for i=1:n
                states(:,1+i)=state+sqrt(lambda+n)*L(:,i);
                states(:,n+i+1)=state-sqrt(lambda+n)*L(:,i);  
            end
            for i=1:2*n+1
                states(:,i)=F*states(:,i) + B*u_p;
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
            % Predict Measurement From Propagated Sigma Points
            measures=zeros(mnum,2*n+1);
            for i=1:2*n+1
                measures(:,i)=height(states(1,i), states(2,i));
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
                    measures(:,i)=height(states(1,i), states(2,i));
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
            n=2;
            mnum=1;
            states=zeros(n,n*2);
            for i=1:n
                states(:,i)=state+sqrt(n)*L(:,i);
                states(:,n+i)=state-sqrt(n)*L(:,i);  
            end
            state=0;
            for i=1:2*n
                states(:,i)=F*states(:,i) + B*u_p;
                state=state+states(:,i)/n/2;
            end
            P=Q;
            for i=1:2*n
                P=P+(state-states(:,i))*(state-states(:,i)).'/n/2;
            end
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
                measures(:,i) = height(states(1,i), states(2,i)); % predicted measurement
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
                    measures(:,i) = height(states(1,i), states(2,i)); % predicted measurement
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
            % predicted state estimate
            x_est(:,k) = F*x_est(:,k-1) + B*u_p;
            n=length(x_est(:,k));
            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update
            mnum=1;
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=ddheight(x_est(1,k),x_est(2,k));

            % measurement matrix
            H = dheight(x_est(1,k),x_est(2,k));
            S=H*P*H.'+R;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*P*Hxx(:,:,jj)*P));
                end
            end
            % Kalman gain
            K=P*H.'*S^(-1);
            % predicted meausrement
            z_p = height(x_est(1,k), x_est(2,k));
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
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=ddheight(x_est(1,k),x_est(2,k));
                H2 = dheight(x_est(1,k),x_est(2,k));
                S2=H2*P*H2.'+R;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*P*Hxx2(:,:,jj)*P));
                    end
                end
                temp=P+K*S2*K.'-P*H2.'*K.'-K*H2*P;
                if(sum(diag(temp))<sum(diag(P)))
                    P=temp;
                else
                    x_est(:,k)=x_est0;
                end
            else
                P = P-K*S*K.';
            end
        else
            %2-IEKF
            % predicted state estimate
            x_est(:,k) = F*x_est(:,k-1) + B*u_p;
            n=length(x_est(:,k));
            % predicted error covariance
            P = F*P*F' + Q;

            %%% Update
            mnum=1;
            Hxx=zeros(n,n,mnum);
            Hxx(:,:,1)=ddheight(x_est(1,k),x_est(2,k));

            % measurement matrix
            H = dheight(x_est(1,k),x_est(2,k));
            S=H*P*H.'+R;
            for ii=1:mnum
                for jj=1:mnum
                    S(ii,jj)=S(ii,jj)+0.5*sum(diag(Hxx(:,:,ii)*P*Hxx(:,:,jj)*P));
                end
            end
            % Kalman gain
            K=P*H.'*S^(-1);
            % predicted meausrement
            z_p = height(x_est(1,k), x_est(2,k));
            % measurement residual
            residual=z-z_p;
            for ii=1:mnum
                residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*P));
            end
            % updated state estimate
            x_est0 = x_est(:,k);
            x_est(:,k) = x_est(:,k) + K*residual;
            %IEKF
            steplennow = norm(K*residual);
            change = 1;
            iter = 1;
            while(change>0.001 && iter < 1000)
                iter = iter + 1;
                z_p = height(x_est(1,k), x_est(2,k));
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=ddheight(x_est(1,k),x_est(2,k));
                residual = z - z_p;
                for ii=1:mnum
                    residual(ii)=residual(ii)-0.5*(x_est0-x_est(:,k)).'*Hxx2(:,:,ii)*(x_est0-x_est(:,k));
                    residual(ii)=residual(ii)-0.5*sum(diag(Hxx(:,:,ii)*P));
                end
                H2 = dheight(x_est(1,k),x_est(2,k));
                S2=H2*P*H2.'+R;
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
                else
                    change = max(abs(dx./x_est(:,k)));
                    x_est(:,k) = x_est(:,k) + dx;
                    K = K2;
                    S = S2;
                end
            end
            % updated error covariance
            if(improvement==1)
                Hxx2=zeros(n,n,mnum);
                Hxx2(:,:,1)=ddheight(x_est(1,k),x_est(2,k));
                H2 = dheight(x_est(1,k),x_est(2,k));
                S2=H2*P*H2.'+R;
                for ii=1:mnum
                    for jj=1:mnum
                        S2(ii,jj)=S2(ii,jj)+0.5*sum(diag(Hxx2(:,:,ii)*P*Hxx2(:,:,jj)*P));
                    end
                end
                temp=P+K*S2*K.'-P*H2.'*K.'-K*H2*P;
                if(sum(diag(temp))<sum(diag(P)))
                    P=temp;
                else
                    x_est(:,k)=x_est0;
                end
            else
                P = P-K*S*K.';
            end
        end
        temp = toc;
        Runtime=Runtime + temp/M;
        P_diag(:,k) = diag(P);
        
        res_d_err(k,m) = norm(x_est(:,k) - x_true(:,k));
    end
    
    res_x_est(:,:,m) = x_est;
    res_x_err(:,:,m) = x_est - x_true;
end

%% get result statistics

x_est_avg = mean(res_x_est,3);
x_err_avg = mean(res_x_err,3);
x_RMSE = zeros(2,N+1); % root mean square error
for k = 1:1:N+1
    x_RMSE(1,k) = sqrt(mean(res_x_err(1,k,:).^2,3));
    x_RMSE(2,k) = sqrt(mean(res_x_err(2,k,:).^2,3));
end
time = (0:1:N)*dt;
end