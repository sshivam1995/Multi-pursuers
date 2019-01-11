function zdot=dzdt(z,u,V)
    
    %beta=atan(lr/(lf+lr)*tan(u(1)));
    
    zdot(1,1)= V*cos(z(3));
    zdot(2,1)= V*sin(z(3));
    zdot(3,1)= u;
end