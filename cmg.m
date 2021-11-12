%% Singularity Avoidance of CMG in Attitude Control of Satellite
clear all
close all
clc

%% Initial Value
w(:,1) = deg2rad([ 0 0 0 ])'; %deg to rad/s (max= 0.1rad/s)
euler_0 = deg2rad( [ 0 0 0 ] ); %deg to rad
euler_cmd = deg2rad( [ -36 20 30 ] ); % deg to rad
q_cur(:,1) = angle2quat(euler_0(1),euler_0(2),euler_0(3), 'ZYX'); % euler to quaternion
q_cmd(:,1) = angle2quat(euler_cmd(1),euler_cmd(2),euler_cmd(3), 'ZYX'); % euler to quaternion
Isc = [0.028 0 0 ; 0 0.028 0 ; 0 0 0.034]; %kgm^2
beta=deg2rad(53.13); %deg to rad
h_a= 1.312*10^(-6)*(100*pi); %Nmsec
zeta = 0.707; %damping ratio
gamma(:,1) = deg2rad([ 0 0 0 0 ]); %deg to rad
ts= 15; % settling time sec
wn=4/(zeta*ts); %natural frequency
k_D=2*zeta*wn*Isc; %D gain
k_P=2*wn^2*Isc; %P gain
dt=0.1; %sec
time=0:dt:50;
limit_GV=pi ; %rad/s
tt=0;
e_0=0.01;

%%
for t=1:length(time)-1
    %% Command Torque
    q_c=[q_cmd(4) q_cmd(3) -q_cmd(2) -q_cmd(1); % skew vector of command quaternion
        -q_cmd(3) q_cmd(4) q_cmd(1) -q_cmd(2);
        q_cmd(2) -q_cmd(1) q_cmd(4) -q_cmd(3);
        q_cmd(1) q_cmd(2) q_cmd(3) q_cmd(4)];

    q_e(:,t)=q_c*q_cur(:,t); %quaternion error
    q_e3(:,t)=[q_e(1,t) q_e(2,t) q_e(3,t)]';
    u_cmd(:,t) = - ( k_P*q_e3(:,t) + k_D*w(:,t) );

    %% Angular Momentum H
    A = [-cos(beta)*sin(gamma(1,t))   -cos(gamma(2,t))             cos(beta)*sin(gamma(3,t))      cos(gamma(4,t));...
      cos(gamma(1,t))            -cos(beta)*sin(gamma(2,t))   -cos(gamma(3,t))               cos(beta)*sin(gamma(4,t));...
      sin(beta)*sin(gamma(1,t)) sin(beta)*sin(gamma(2,t))   sin(beta)*sin(gamma(3,t))    sin(beta)*sin(gamma(4,t)) ];

    h_cmgxyz(:,t) = [A(1,1)+A(1,2)+A(1,3)+A(1,4);...
             A(2,1)+A(2,2)+A(2,3)+A(2,4);...
             A(3,1)+A(3,2)+A(3,3)+A(3,4)]*h_a;

    %% H_dot  
    u_cmgxyz(:,t)=-( u_cmd(:,t) + cross(w(:,t),h_cmgxyz(:,t)) );
   
    %% Jacobian Matrix
    dotA = [ -cos(beta)*cos(gamma(1,t)) sin(gamma(2,t)) cos(beta)*cos(gamma(3,t)) -sin(gamma(4,t));...
           -sin(gamma(1,t)) -cos(beta)*cos(gamma(2,t)) sin(gamma(3,t)) cos(beta)*cos(gamma(4,t));...
           sin(beta)*cos(gamma(1,t)) sin(beta)*cos(gamma(2,t)) sin(beta)*cos(gamma(3,t)) sin(beta)*cos(gamma(4,t))];
       
    %% Gimbal velocity
    tt=tt+dt;
    m(:,t) = det(dotA*dotA');
    a = e_0*exp(-10*m(:,t));
    E=[          1             e_0*sin(tt+pi)    e_0*sin(tt+pi/2);...
         e_0*sin(tt+pi)            1             e_0*sin(tt);...
         e_0*sin(tt+pi/2)      e_0*sin(tt)             1           ];
    dotAA = dotA'*inv( dotA*dotA' + a*E );
    dotgamma(:,t) = dotAA*u_cmgxyz(:,t)/h_a;
    
    for i=1:4
        if dotgamma(i,t) > limit_GV
            dotgamma_R(i,t)=limit_GV;
        elseif dotgamma(i,t) < -limit_GV
            dotgamma_R(i,t)=-limit_GV;
        else
            dotgamma_R(i,t)=dotgamma(i,t);
        end
    end
    
    gamma(:,t+1)=gamma(:,t)+dotgamma_R(:,t)*dt;
    
    %% Real torque
    u_real(:,t) = -( h_a*dotA*dotgamma_R(:,t) + cross(w(:,t),h_cmgxyz(:,t)));
       
    %% Dynamic Equation
    u_real_ode(:,1) = u_real(:,t);
    w_ode(:,1) = w(:,t);   
    [T1(:,t),W]=ode45(@(t1,w_ode) Isc\( u_real_ode(:,1)  -cross(w_ode(:,1),Isc*w_ode(:,1)) ), [time(t):dt/2:time(t+1)] ,[w_ode(:,1)]);
    w(:,t+1) = W(3,:);
 
    %% Kinematics Equation
    q_w=[0 w(3,t) -w(2,t) w(1,t); -w(3,t) 0 w(1,t) w(2,t);...
    w(2,t) -w(1,t) 0 w(3,t); -w(1,t) -w(2,t) -w(3,t) 0 ];
    q_cur_ode(:,1) = q_cur(:,t);
    [T2(:,t),Q]=ode45(@(t,q_cur_ode) 0.5*q_w*q_cur_ode(:,1),[time(t):dt/2:time(t+1)],[q_cur_ode(:,1)]);
    q_cur(:,t+1) = Q(3,:);
    
    %% Quaternion -> euler angle(deg)
    RPY(:,t) = rad2deg(quat2eul( [q_cur(1,t) q_cur(2,t) q_cur(3,t) q_cur(4,t) ],'ZYX' )); 
  
end


%% Plot
% figure(1)
% plot(T1(1,:), RPY)
% title('Spacecraft Euler Angle','fontsize',15);
% legend({'Roll','Pitch','Yaw'},'fontsize',15)
% xlabel('sec','fontsize',15);
% ylabel('degree','fontsize',15);
% grid on;
% hold on;

% figure(2)
% plot(time, w)
% legend({'\omega_x','\omega_y','\omega_z'},'fontsize',15)
% title('Spacecraft Angular Velocity','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('rad/s','fontsize',15);
% grid on;
% hold on;
% 
% figure(3)
% plot(T1(1,:), dotgamma_R)
% legend({'cmg_c_l_u_s_t_e_r_1','cmg_c_l_u_s_t_e_r_2','cmg_c_l_u_s_t_e_r_3','cmg_c_l_u_s_t_e_r_4'},'fontsize',15)
% title('Gimbal velocity','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('rad/s','fontsize',15);
% grid on;
% hold on;
% 
% figure(4)
% plot(time, gamma*180/pi)
% legend({'cmg_c_l_u_s_t_e_r_1','cmg_c_l_u_s_t_e_r_2','cmg_c_l_u_s_t_e_r_3','cmg_c_l_u_s_t_e_r_4'},'fontsize',15)
% title('Gimbal angle \gamma','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('deg','fontsize',15);
% grid on;
% hold on;
% 
% figure(5)
% plot( T1(1,:), u_real )
% legend({'\tau_x', '\tau_y', '\tau_z'},'fontsize',15)
% title('CMG Torque','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('Nm','fontsize',15);
% grid on;
% hold on;
% 
% figure(6)
% plot( T1(1,:), m )
% title('Singularity State','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('det(A^\primeA^\prime^T)','fontsize',15);
% grid on;
% hold on;

%% Comparison plot  // ¸ÅÆ®·¦ v.s. ÆÄÀÌ½ã
% RPY_Py = xlsread('RPY.csv') 
% 
% figure(1)
% plot(T1(1,:), RPY_Py(:,1), '--b', T1(1,:), RPY_Py(:,2), '--r',T1(1,:), RPY_Py(:,3), '--y','LineWidth',1)
% legend({'Roll_M_A_T_L_A_B','Pitch_M_A_T_L_A_B','Yaw_M_A_T_L_A_B',...
%         'Roll_P_y_t_h_o_n','Pitch_P_y_t_h_o_n','Yaw_P_y_t_h_o_n'},'fontsize',15)
% grid on;

%% Comparison plot // ÆÄÀÌ½ã v.s. ÆÄÀÌ½ã_ÇÇµå¹é
% RPY_Py = xlsread('RPY_.csv');
% RPY_Py_close = xlsread('RPY_close.csv');
% w_Py = xlsread('w_.csv');
% w_Py_close = xlsread('w_close.csv');
% dotgamma_real_Py = xlsread('dotgamma_real_.csv');
% dotgamma_real_Py_close = xlsread('dotgamma_real_close.csv');
% gamma_Py = xlsread('gamma_.csv');
% gamma_Py_close = xlsread('gamma_close.csv');
% u_real_Py = xlsread('u_real_.csv');
% u_real_Py_close = xlsread('u_real_close.csv');
% m_Py = xlsread('m_.csv');
% m_Py_close = xlsread('m_close.csv');
% 

figure(3)
plot(time, dotgamma_real_Py(:,1))
legend({'cmg_c_l_u_s_t_e_r_1'},'fontsize',11)
title('Gimbal velocity','fontsize',15);
xlabel('sec','fontsize',15);
ylabel('rad/s','fontsize',15);
grid on;
hold on;

figure(3)
plot(time, dotgamma_real_Py_close(:,1), 'r','LineWidth',1)
legend({'cmg_c_l_u_s_t_e_r_1','cmg_c_l_u_s_t_e_r_1_,_c_l_o_s_e'},'fontsize',15)
grid on;

% 
% figure(1)
% plot(time, RPY_Py)
% title('Spacecraft Euler Angle','fontsize',15);
% legend({'Roll','Pitch','Yaw'},'fontsize',15)
% xlabel('sec','fontsize',15);
% ylabel('degree','fontsize',15);
% grid on;
% hold on;
% 
% figure(1)
% plot(time, RPY_Py_close(:,1), '--b', time, RPY_Py_close(:,2), '--r',time, RPY_Py_close(:,3), '--y','LineWidth',1)
% legend({'Roll','Pitch','Yaw',...
%         'Roll_c_l_o_s_e','Pitch_c_l_o_s_e','Yaw_c_l_o_s_e'},'fontsize',15)
% grid on;
% 
% figure(2)
% plot(time, w_Py)
% legend({'\omega_x','\omega_y','\omega_z'},'fontsize',15)
% title('Spacecraft Angular Velocity','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('rad/s','fontsize',15);
% grid on;
% hold on;
% 
% figure(2)
% plot(time, w_Py_close(:,1), '--b', time, w_Py_close(:,2), '--r',time, w_Py_close(:,3), '--y','LineWidth',1)
% legend({'\omega_x','\omega_y','\omega_z',...
%         '\omega_x_,_c_l_o_s_e','\omega_y_,_c_l_o_s_e','\omega_z_,_c_l_o_s_e'},'fontsize',15)
% grid on;
% 
% figure(3)
% plot(time, dotgamma_real_Py)
% legend({'cmg_c_l_u_s_t_e_r_1','cmg_c_l_u_s_t_e_r_2','cmg_c_l_u_s_t_e_r_3','cmg_c_l_u_s_t_e_r_4'},'fontsize',11)
% title('Gimbal velocity','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('rad/s','fontsize',15);
% grid on;
% hold on;
% 
% figure(3)
% plot(time, dotgamma_real_Py_close(:,1), '--b', time, dotgamma_real_Py_close(:,2), '--r',...
%       time, dotgamma_real_Py_close(:,3), '--y',time, dotgamma_real_Py_close(:,4), '--m','LineWidth',1)
% legend({'cmg_c_l_u_s_t_e_r_1','cmg_c_l_u_s_t_e_r_2','cmg_c_l_u_s_t_e_r_3','cmg_c_l_u_s_t_e_r_4',...
%         'cmg_c_l_u_s_t_e_r_1_,_c_l_o_s_e','cmg_c_l_u_s_t_e_r_2_,_c_l_o_s_e','cmg_c_l_u_s_t_e_r_3_,_c_l_o_s_e','cmg_c_l_u_s_t_e_r_4_,_c_l_o_s_e'},'fontsize',11)
% grid on;
% 
% figure(4)
% plot(time, gamma_Py*180/pi)
% legend({'cmg_c_l_u_s_t_e_r_1','cmg_c_l_u_s_t_e_r_2','cmg_c_l_u_s_t_e_r_3','cmg_c_l_u_s_t_e_r_4'},'fontsize',11)
% title('Gimbal angle \gamma','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('deg','fontsize',15);
% grid on;
% hold on;
% 
% figure(4)
% plot(time, gamma_Py_close(:,1)*180/pi, '--b', time, gamma_Py_close(:,2)*180/pi, '--r',...
%       time, gamma_Py_close(:,3)*180/pi, '--y',time, gamma_Py_close(:,4)*180/pi, '--m','LineWidth',1)
% legend({'cmg_c_l_u_s_t_e_r_1','cmg_c_l_u_s_t_e_r_2','cmg_c_l_u_s_t_e_r_3','cmg_c_l_u_s_t_e_r_4',...
%         'cmg_c_l_u_s_t_e_r_1_,_c_l_o_s_e','cmg_c_l_u_s_t_e_r_2_,_c_l_o_s_e','cmg_c_l_u_s_t_e_r_3_,_c_l_o_s_e','cmg_c_l_u_s_t_e_r_4_,_c_l_o_s_e'},'fontsize',11)
% grid on;
% 
% figure(5)
% plot( time, u_real_Py )
% legend({'\tau_x', '\tau_y', '\tau_z'},'fontsize',15)
% title('CMG Torque','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('Nm','fontsize',15);
% grid on;
% hold on;
% 
% figure(5)
% plot(time, u_real_Py_close(:,1), '--b', time, u_real_Py_close(:,2), '--r',time, u_real_Py_close(:,3), '--y','LineWidth',1)
% legend({'\tau_x', '\tau_y', '\tau_z',...
%         '\tau_x_,_c_l_o_s_e', '\tau_y_,_c_l_o_s_e', '\tau_z_,_c_l_o_s_e'},'fontsize',15)
% grid on;
% 
% figure(6)
% plot( time, m_Py )
% title('Singularity State','fontsize',15);
% xlabel('sec','fontsize',15);
% ylabel('det(A^\primeA^\prime^T)','fontsize',15);
% grid on;
% hold on;
% 
% figure(6)
% plot(time, m_Py_close(:,1), '--r','LineWidth',1)
% legend({'m', 'm_c_l_o_s_e'},'fontsize',15)
% grid on;