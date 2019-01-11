close all, clear all

%% goal defined
goal=[15;0];

%% sim parameters

dt=0.01; tf=400; t=0:dt:tf;
N=size(t,2);
alpha=30;
beta3_max=5;
beta1=1; beta2=1; beta3=beta3_max;
%% evader and pursuer parameters

V_E=0.1;
gamma=1/3;      % ratio of evader to pursuer speed
V_P=V_E/gamma;

% R_c/R_l must be greater than sqrt(1-gamma^2)+asin(gamma)-1 = 0.2826 for capture to be possible. We pick 0.3 

scale=5;
R_t=scale*0.2;            % Turn radius
R_c=R_t*0.2;         % Capture radius

udot=zeros(2,N);
angular_vel_sat=V_P/R_t;

%% States of evader and pursuer
Z_P1=zeros(3,N);
Z_P1(:,1)=[-5;-0;0.1];

Z_P2=zeros(3,N);
Z_P2(:,1)=[-8;-4;0.1];

u_P=zeros(2,N);

Z_E=zeros(2,N);
u_E=u_P;

%% Track vehicle motion
pos_tracker_evader=zeros(2,N);
pos_tracker_evader(:,1)=Z_E(1:2,1);

pos_tracker_pursuer=zeros(4,N);
pos_tracker_pursuer(1:2,1)=Z_P1(1:2,1);
pos_tracker_pursuer(3:4,1)=Z_P2(1:2,1);

distance=zeros(3,N);

%% derivative of z at all time instances
zdot_P1=zeros(3,N); zdot_P2=zeros(3,N); zdot_E=zeros(2,N);

g_u_P1_pursue=zeros(1,N);
g_u_prime_P1=zeros(1,1,N);

g_u_P2_pursue=zeros(1,N);
g_u_prime_P2=zeros(1,1,N);

%% main loop
for i=1:N 
    %% Check capture condition
    distance(1,i)= norm(pos_tracker_pursuer(1:2,i)-pos_tracker_evader(:,i));
    distance(2,i)= norm(pos_tracker_pursuer(3:4,i)-pos_tracker_evader(:,i));
    distance(3,i)= norm(pos_tracker_pursuer(1:2,i)-pos_tracker_pursuer(3:4,i));
%     if distance(i)<R_c
%         Winner=1;
%         break
%     end
    
    %% Check goal reach condition
    if norm(goal-pos_tracker_evader(:,i))<0.001
        Winner=2;
 
    end
    
    %% Set horizon depending on distance
    if distance(1,i)>2*R_t
       T1=0.8;
    else
       T1=0.8;
    end
    delta_T1=0.1*T1;

    if distance(2,i)>2*R_t
       T2=0.5;
    else
       T2=0.3;
    end
    delta_T2=0.1*T2;
    
     %% set weight for go to goal behaviour dependent on distance
     min_distance=min(distance(1,i),distance(2,i));
     max_distance=max(distance(1,i),distance(2,i));
     
     if max_distance>2*R_t
        beta3=0;
        beta2=2;
     else
        beta3=beta3_max-beta3_max*max_distance/(2*R_t);
        beta2=1;
     end    
         
    if min_distance>2*R_t        
        W=1;
        beta1=0.2;
      %  beta3=0;
    elseif min_distance>R_t && min_distance<=2*R_t
        W=(min_distance-R_t)/R_t;
        beta1=0.5;
    else
        W=0;
        beta1=1;
    end
    
    if norm(goal-pos_tracker_evader(:,i))<R_t*gamma
        W=1; 
    end    
    
    %% set weights for circling and pursuit behaviours
    
 %   if distance(1,i)<=2*R_t
    
     %% Future prediction
 
    if distance(1,i)<distance(2,i)
        closest_pursuer=1;
    else    
        closest_pursuer=2;
    end    
    
     homing_vector=pos_tracker_evader(:,i)-pos_tracker_pursuer(2*closest_pursuer-1:2*closest_pursuer,i);
     homing_direction=atan2((homing_vector(2)),(homing_vector(1)));
     
    if min_distance>R_t
        vel_angle_evader_optimal=homing_direction;
        mode=1;
    else
%% decide evasion depending on closest pursuer, we have 2 choices, angle+pi/2 or angle-pi/2; 
% we choose the one that gets evader closer to goal; once choice is made,
% we fix that till we're out of closest distance<= turning radius
        
        if closest_pursuer==1
            if mode==2
                vel_angle_evader_optimal=mod(Z_P1(3,i),2*pi)+pi/2;
            elseif mode==3
                vel_angle_evader_optimal=mod(Z_P1(3,i),2*pi)-pi/2;
            else
                angle_goal=atan2(goal(2)-pos_tracker_pursuer(2,i),goal(1)-pos_tracker_pursuer(1,i));
                angle_pursuer=mod(Z_P1(3,i),2*pi);
                
                if angle_pursuer>pi
                   angle_pursuer=angle_pursuer-2*pi;
                end
                
                if angle_goal>angle_pursuer
                    vel_angle_evader_optimal=mod(Z_P1(3,i),2*pi)+pi/2;
                    mode=2;
                else
                    vel_angle_evader_optimal=mod(Z_P1(3,i),2*pi)-pi/2;
                    mode=3;
                end    
            end
            
        else
            if mode==2
                vel_angle_evader_optimal=mod(Z_P2(3,i),2*pi)+pi/2;
            elseif mode==3
                vel_angle_evader_optimal=mod(Z_P2(3,i),2*pi)-pi/2;
            else
                angle_goal=atan2(goal(2)-pos_tracker_pursuer(4,i),goal(1)-pos_tracker_pursuer(2,i));
                angle_pursuer=mod(Z_P2(3,i),2*pi);
                
                if angle_pursuer>pi
                   angle_pursuer=angle_pursuer-2*pi;
                end
                
                if angle_goal>angle_pursuer
                    vel_angle_evader_optimal=mod(Z_P2(3,i),2*pi)+pi/2;
                    mode=2;
                else
                    vel_angle_evader_optimal=mod(Z_P2(3,i),2*pi)-pi/2;
                    mode=3;
                end    
            end
        end
    end    

%% for boundary conditions    
    if i>1
        old_input_pursuer1=u_P(1,i-1);
        old_input_pursuer2=u_P(2,i-1);
        
        old_velocity1= zdot_P1(1:2,i-1);
        old_velocity2= zdot_P2(1:2,i-1);
    else
        old_input_pursuer1=0;
        old_input_pursuer2=0;
        
        old_velocity1= 0;
        old_velocity2= 0;
    end
    
 %   vel_angle_evader_optimal=0;
    
%% calculate g and dgdu    
    [g_u,g_u_prime_P1,g_u_prime_P2]=g_rt_integral(Z_P1(:,i),old_input_pursuer1, Z_P2(:,i),old_input_pursuer2, T1, delta_T1,V_P, pos_tracker_evader(:,i), vel_angle_evader_optimal, V_E, beta1,beta2,beta3,R_t);
 %   [g_u_P2_pursue,g_u_prime_P2_pursue]=g_rt2(Z_P2(:,i),old_input_pursuer2, T2, delta_T2,V_P, pos_tracker_evader(:,i), vel_angle_evader_optimal, V_E);

    pursuit_ang_vel1=old_input_pursuer1+(1/(g_u_prime_P1))*(0-g_u)*alpha*dt;
    pursuit_ang_vel2=old_input_pursuer2+(1/(g_u_prime_P2))*(0-g_u)*alpha*dt;
    
%     [g_u_P1_circle,g_u_prime_P1_circle] = g_rt3(Z_P1(:,i),old_input_pursuer1, T1, delta_T1,V_P, pos_tracker_pursuer(3:4,i), old_velocity2);
%     [g_u_P2_circle,g_u_prime_P2_circle] = g_rt3(Z_P2(:,i),old_input_pursuer2, T2, delta_T2,V_P, pos_tracker_pursuer(1:2,i), old_velocity1);
% 
%     circle_ang_vel1=old_input_pursuer1+(1/(g_u_prime_P1_circle))*(2*R_t-g_u_P1_circle)*alpha*dt;
%     circle_ang_vel2=old_input_pursuer2+(1/(g_u_prime_P2_circle))*(2*R_t-g_u_P2_circle)*alpha*dt;  

    
    u_P(1,i)= pursuit_ang_vel1;
    u_P(2,i)= pursuit_ang_vel2;
%% apply input saturation due to turning radius
    if abs(u_P(1,i))>angular_vel_sat
        u_P(1,i)=u_P(1,i)*angular_vel_sat/abs(u_P(1,i));
    end   
    
    if abs(u_P(2,i))>angular_vel_sat
        u_P(2,i)=u_P(2,i)*angular_vel_sat/abs(u_P(2,i));
    end
    
%% Assign velocity to evader and pursuer and update states
    
    zdot_P1(:,i)=dzdt(Z_P1(:,i),u_P(1,i),V_P);  
    zdot_P2(:,i)=dzdt(Z_P2(:,i),u_P(2,i),V_P);
    
    Z_P1(:,i+1)=Z_P1(:,i)+zdot_P1(:,i)*dt;  
    pos_tracker_pursuer(1:2,i+1)=[Z_P1(1,i+1);Z_P1(2,i+1)];
    
    Z_P2(:,i+1)=Z_P2(:,i)+zdot_P2(:,i)*dt;  
    pos_tracker_pursuer(3:4,i+1)=[Z_P2(1,i+1);Z_P2(2,i+1)]; 
    
    vel_angle_evasion=normrnd(vel_angle_evader_optimal,0.0);     % normal distribution of evasion angle around optimal value
    vel_evasion=[V_E*cos(vel_angle_evasion);
                 V_E*sin(vel_angle_evasion)];
    
    vel_gtg=(goal-pos_tracker_evader(:,i))*V_E/norm(goal-pos_tracker_evader(:,i));
             
    zdot_E_direction=(1-W)*vel_evasion+W*vel_gtg;
    
    zdot_E(:,i)=zdot_E_direction*V_E/norm(zdot_E_direction);
    
    Z_E(:,i+1)=Z_E(:,i)+zdot_E(:,i)*dt;
    pos_tracker_evader(:,i+1)=[Z_E(1,i+1);Z_E(2,i+1)];
        
end    

% for i=1:N-1
%     vel_track(1,i)=norm(path(:,i)-path(:,i+1))/dt;
% end
% figure (1)
% plot (t(1,1:N-1),vel_track)

% Plot of accln    
% figure (1)
% plot(t,u(2,:))

% 
% 
% % Plot of jerk
% figure (2)
% plot(t,udot(2,:))
% title('Longitudinal jerk vs time');
% 

% PLot of tracking error norm
% figure (3)
% 
% error=pos_tracker-r;   
% 
% norm_error=vecnorm(error);
% 
% plot (t,norm_error);
% title('tracking error norm vs time');
% 
% 
% figure (4)
% 
% path_heading=mod(atan2(V(2,1:N),V(1,1:N)),2*pi);
% actual_heading=mod(Z(6,1:N),2*pi);
% 
% error_heading=mod(actual_heading-path_heading,2*pi);
% for i=1:size(error_heading,2)
%         if error_heading(i)>6
%             error_heading(i)=error_heading(i)-2*pi;
%         end 
% end
% 
% plot(t,error_heading)
% title('Heading error vs time');
% 
% 
% 
% figure (5)
% 
% [close_point,lateral_error_norm,arc]=distance2curve(r',pos_tracker');
% 
% plot(t,lateral_error_norm)
% title ('Normal Error')
%  

figure (1)

plot(0:dt:(i-1)*dt,distance(1,1:i))
hold on
plot(0:dt:(i-1)*dt,distance(2,1:i))
cut_off=R_c*ones(1,i);
%plot(0:dt:(i-1)*dt,cut_off)
radius_turn=R_t*ones(1,i);
plot(0:dt:(i-1)*dt,radius_turn)
title('Distance of pursuers from evader')
legend('Pursuer 1', 'Pursuer 2', 'turning radius')


figure (6)
plot(pos_tracker_pursuer(1,1:i),pos_tracker_pursuer(2,1:i), '--r');
hold on
plot(pos_tracker_pursuer(3,1:i),pos_tracker_pursuer(4,1:i), '--y');
 plot (pos_tracker_evader(1,1:i),pos_tracker_evader(2,1:i),'b');
% title ('Tracking path')
% % xlim([-20 100]) 
% % ylim([-20 100])


goal_pos=goal*ones(1,i);

title('Trajectories of pursuers and evader')
xlim([-10 16])
ylim([-10 16])
s = plot(pos_tracker_pursuer(1,1),pos_tracker_pursuer(2,1),'o','MarkerFaceColor','red');      % bot 1
t = plot(pos_tracker_pursuer(3,1),pos_tracker_pursuer(4,1),'o','MarkerFaceColor','yellow');
q = plot(pos_tracker_evader(1,1),pos_tracker_evader(2,1),'o','MarkerFaceColor','blue');                          % ref
r = plot(goal_pos(1,1),goal_pos(2,1),'o','MarkerFaceColor','black');

legend('Pursuer 1 trajectory','Pursuer 2 trajectory','Evader trajectory', 'Pursuer 1', 'Pursuer 2', 'Evader', 'goal');
for k = 11:10:i-mod(i,10)+1
    s.XData = pos_tracker_pursuer(1,k);
    s.YData = pos_tracker_pursuer(2,k);
    
    t.XData = pos_tracker_pursuer(3,k);
    t.YData = pos_tracker_pursuer(4,k);
%       
     q.XData = pos_tracker_evader(1,k);
     q.YData = pos_tracker_evader(2,k);
%     
    drawnow
end
    

% 
% figure (7)
% plot(t,r(2,:))
% hold on
% plot(t,Z(5,:))
% legend('reference','actual');
% title ('Y position')



