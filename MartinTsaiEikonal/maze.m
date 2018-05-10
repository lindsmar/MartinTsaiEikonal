function [F]=maze(N,K)
% N=10;
% K=50;
M=N*K+1;
F=ones(M,M);
h=1/(N*K);

x1=.4;
y1=.18;
x2=.1;
y2=.605;
x3=.75;
y3=.45;
x4=.75;
y4=.15;
x5=.4;
y5=.82;
x6=.7;
y6=.395;
for i=1:M
    for j=1:M
        x=(j-1)*h;
        y=(N*K+1-i)*h;
        R=sqrt((x-x1)^2+(y-y1)^2);
        theta=atan2(y-y1,x-x1);

        
        if (theta>=2*pi/2.9 || theta<=-pi+pi/2.9)
            if R>=.25 && R<=.27
                F(i,j)=1000;
            end
        end
        R=sqrt((x-x2)^2+(y-y2)^2);
        theta=atan2(y-y2,x-x2);

        
        if (theta>=-pi+2*pi/2.9 && theta<=pi/2.9)
            if R>=.25 && R<=.27
                F(i,j)=1000;
            end
        end
        
        R=sqrt((x-x5)^2+(y-y5)^2);
        theta=atan2(y-y5,x-x5);

        
        if (theta>=-pi+2*pi/2.9 && theta<=pi/2.9)
            if R>=.25 && R<=.27
                F(i,j)=1000;
            end
        end
        
         R=sqrt((x-x6)^2+(y-y6)^2);
        theta=atan2(y-y6,x-x6);

        
        if (theta>=2*pi/2.9 || theta<=-pi+pi/2.9)
            if R>=.25 && R<=.27
                F(i,j)=1000;
            end
        end
        
        
        R=sqrt((x-x3)^2+(y-y3)^2);
       
        if R<=.08
            F(i,j)=1000;
        end
        
        R=sqrt((x-x4)^2+(y-y4)^2);
       
        if R<=.04
            F(i,j)=.01;
        end
        
      
        
        
                

        
    end
end
        

end