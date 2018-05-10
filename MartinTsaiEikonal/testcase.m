function [ORIG,bdrycond,R]=testcase(N,K,test)

R=ones(51,51);
H=1/N; 
h=1/(N*K);


switch test 
     case 1 %R=1+.99sin(2pix)sin(2piy)
        Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        % source point
        z=1;
        pt=[0.5,0.5];
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        % slowness function
        [X1,X2]=ndgrid(0:h:1,0:h:1);
%                 X1=flipud(X1);
%         X2=flipud(X2);

        R=1 + .99*sin(2*pi*X1).*sin(2*pi*X2);

        % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,50,h);
        
    case 2 %R=1+.5sin(20pix)sin(20piy)
    
        Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        [X1,X2]=ndgrid(0:h:1,0:h:1);
%         X1=flipud(X1);
%         X2=flipud(X2);
        % source point
        z=1;
        pt=[0.5,0.5];
         
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        % slowness function
        R=1 + .5*sin(50*pi*X1).*sin(50*pi*X2);

        % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,50,h);
        
    case 3
        Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        % source point
        z=1;
        pt=[0.5,0.5];
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        % slowness function
         R=ones(N*K+1,N*K+1);
        % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,50,h);
    case 4
        Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        z=1;
        pt=[0,0];
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        % slowness function
        R=maze(N,K);

         % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,50,h);
        
    case 5
        Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        z=1;
        pt=[0,0];
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        
        % slowness function
        R=ones(N*K+1,N*K+1);
        x1min=ceil(.26/h)+1;
        x1max=ceil(.27/h)+1;
        ytop=N*K+1-ceil(.6/h);
        R(ytop:end,x1min:x1max)=0.01;
        % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,50,h);
    
    case 6 %generalization of case 2
         Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        [X1,X2]=ndgrid(0:h:1,0:h:1);
        X1=flipud(X1);
        X2=flipud(X2);
        % source point
%         z=2;
%         pt=[0.25,0.25; .75,.75];
       z=2;
       pt=[0.35,0.35;0.65,0.65];
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        % slowness function
        epsilon=(abs(X1)+abs(X2)+.001)/50;
        R=1 + .5*sin(pi*X1./epsilon).*sin(pi*X2./epsilon);

        % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,50,h);
    case 7 %homogenization
        Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        
        % source point
        z=1;
        pt=[0.5,0.5];
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        % slowness function 
        R=makespeed(7*h,N*K+1,2,3);
        % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,100,h);
    case 8 %random coefficients
        Uold=1000*ones(N*K+1,N*K+1);
        BdryCond=-1*ones(N*K+1,N*K+1);
        obs=zeros(N*K+1,N*K+1);
        % source point
        z=1;
        pt=[0.5,0.5];
        for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
            BdryCond(m,k)=0;
        end
        % slowness function
        R=makespeed(7*h,N*K+1,2,2);
        % overall fine solution
        [ORIG, ~, ~]=fastsweeporig(Uold,R, BdryCond, obs,50,h);
        
end



 
    bdrycond=-1*ones(N*K+1,N*K+1);
    % Set up boundary conditions for new method
    for n=1:z
            k=ceil(pt(n,1)/h)+1;
            m=N*K+1-ceil(pt(n,2)/h);
        if rem(k,K)==1
            if k==1
                a=1;
                s=1;
                c=K+1;
                t=K+1;
            elseif k==N*K+1
                a=N*K+1-K;
                s=a;
                c=N*K+1;
                t=c;
            else
                a=k-K;
                s=a;
                c=k+K;
                t=c;
            end
        else
            a=floor(k/K)*K+1-K;
            s=floor(k/K)*K+1;
            t=floor(k/K)*K+1+K;
            if a<1
                a=1;
            end

            if a+K*3>N*K+1
                c=N*K+1;
            else
                c=a+K*3;
            end
        end

        if rem(m,K)==1
            if m==1
                b=1;
                u=1;
                d=K+1;
                v=K+1;
            elseif m==N*K+1
                b=N*K+1-K;
                u=b;
                d=N*K+1;
                v=d;
            else
                b=m-K;
                u=b;
                d=m+K;
                v=d;
            end
        else
            b=floor(m/K)*K+1-K;
            u=floor(m/K)*K+1;
            v=floor(m/K)*K+1+K;
            if b<1
                b=1;
            end

            if b+K*3>N*K+1
                d=N*K+1;
            else
                d=b+K*3;
            end
        end


            bdrycond(b:d,s:t)=ORIG(b:d,s:t);
            bdrycond(u:v,a:c)=ORIG(u:v,a:c);

    end


end
