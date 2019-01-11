function  [g_u_sum,g_u_prime_1_sum,g_u_prime_2_sum]= g_rt_integral(z0_1,u1,z0_2,u2,T,delta_t,V, initial_pos_evader, evasion_angle, evasion_speed, beta1, beta2,beta3,turning_radius)

M=fix(T/delta_t);
g=4;

z1=zeros(3,M);
z1(:,1)=z0_1;

z2=zeros(3,M);
z2(:,1)=z0_2;

distance=zeros(1,M);
distance(1)=norm(z1(1:2,1)-z2(1:2,1));


velocity=   [evasion_speed*cos(evasion_angle);
            evasion_speed*sin(evasion_angle)];

estimate_position=zeros(2,M);
estimate_position(:,1)=initial_pos_evader+velocity*T/M;

dzdu1=zeros(3,1,M);

dzdu2=zeros(3,1,M);

e1=zeros(1,M);
e2=zeros(1,M);

e1(1)=sqrt((z1(1,1)-estimate_position(1,1))^2+(z1(2,1)-estimate_position(2,1))^2);
e2(1)=sqrt((z2(1,1)-estimate_position(1,1))^2+(z2(2,1)-estimate_position(2,1))^2);

g_u=zeros(1,M);
g_u(1)=beta1*(e1(1)^2*e2(1)^2)/(e1(1)^2+e2(1)^2)+beta2*(e1(1)^2+e2(1)^2);

g_u_sum=g_u(1);

g_u_prime_1_sum=0;
g_u_prime_2_sum=0;

for i=2:M
    
    estimate_position(:,i)=initial_pos_evader+velocity*i*delta_t;
    
    zdot1=dzdt(z1(:,i-1),u1,V);
    z1(:,i)=z1(:,i-1)+zdot1*delta_t;
    
    zdot2=dzdt(z2(:,i-1),u2,V);
    z2(:,i)=z2(:,i-1)+zdot2*delta_t;
    
    distance(i)=norm(z1(1:2,i)-z2(1:2,i));
    
    dfdz1=[0 0 -V*sin(z1(3,i));
           0 0 V*cos(z1(3,i));
           0 0 0];
    dfdu1=[0;0;1];
    
    dfdz2=[0 0 -V*sin(z2(3,i));
           0 0 V*cos(z2(3,i));
           0 0 0];
    dfdu2=[0;0;1];
    
    dzdu1(:,:,i)=dzdu1(:,:,i-1)+(dfdz1*dzdu1(:,:,i-1)+dfdu1)*delta_t;     
    dzdu2(:,:,i)=dzdu2(:,:,i-1)+(dfdz2*dzdu2(:,:,i-1)+dfdu2)*delta_t;     
    
    e1(i)=sqrt((z1(1,i)-estimate_position(1,i))^2+(z1(2,i)-estimate_position(2,i))^2);
    e2(i)=sqrt((z2(1,i)-estimate_position(1,i))^2+(z2(2,i)-estimate_position(2,i))^2);
    
    g_u(i)=beta1*(e1(i)^2*e2(i)^2)/(e1(i)^2+e2(i)^2)+beta2*(e1(i)^2+e2(i)^2)+beta3*exp(-g*distance(i));
    
    
    
    dgde1(i)=2*beta1*(e2(i)^4*e1(i))/(e1(i)^2+e2(i)^2)^2+beta2*2*e1(i);
    dgde2(i)=2*beta1*(e1(i)^4*e2(i))/(e1(i)^2+e2(i)^2)^2+beta2*2*e2(i);
    
    de1du1(i)=(1/e1(i))*((z1(1,i)-estimate_position(1,i))*dzdu1(1,1,i)+(z1(2,i)-estimate_position(2,i))*dzdu1(2,1,i));
    de2du2(i)=(1/e2(i))*((z2(1,i)-estimate_position(1,i))*dzdu2(1,1,i)+(z2(2,i)-estimate_position(2,i))*dzdu2(2,1,i));
    
    d_distancedu1(i)=(1/distance(i))*(z1(1,i)-z2(1,i))*dzdu1(1,1,i)+(z1(2,i)-z2(1,i))*dzdu1(2,1,i);
    d_distancedu2(i)=(1/distance(i))*(z2(1,i)-z1(1,i))*dzdu2(1,1,i)+(z2(2,i)-z1(1,i))*dzdu2(2,1,i);
    
    g_u_prime_1(i)=dgde1(i)*de1du1(i)-g*beta3*exp(-g*distance(i))*d_distancedu1(i);
    g_u_prime_2(i)=dgde2(i)*de2du2(i)-g*beta3*exp(-g*distance(i))*d_distancedu2(i);

    g_u_sum=g_u_sum+g_u(i)*T/M;
    
    g_u_prime_1_sum=g_u_prime_1_sum+g_u_prime_1(i);
    g_u_prime_2_sum=g_u_prime_2_sum+g_u_prime_2(i);
end



end