function [ C] = makespeed(epsilon, N,c0, test )
% test=1;
% N=1001;
% c0=2;
% epsilon=1/10;
epsilon=ceil(epsilon/(1/(N-1)));
remainder=rem(N,epsilon);
C=ones(N,N);

switch test
    case 1 %checkerboard
        C2=ones(N,N);
        epsilon=2*epsilon;
        remainder=rem(N,epsilon);

        for i=2:2*epsilon:N-remainder
            for j=2:2*epsilon:N-remainder
                C(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
            end
        end
        for i=epsilon+2:2*epsilon:N-remainder
            for j=epsilon+2:2*epsilon:N-1-remainder
                C(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
            end
        end
         for i=2:2*epsilon:N-remainder
            for j=2+epsilon:2*epsilon:N-remainder
                C2(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
            end
        end
        for i=epsilon+2:2*epsilon:N-remainder
            for j=2:2*epsilon:N-1-remainder
                C2(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
            end
        end
        
        C1=hto2h(C,N);
        C1=flipud(C1);
        C2=hto2h(C2,N);
        C4=C2(2:end,:);
        C2=fliplr(C2(:,2:end));
        C2=flipud(C2);
        C3=flipud(fliplr(C1(1:end-1,2:end)));
        
          C=[C2,C1;C3,C4];
        
    case 2 %random checkerboard
        
%          C2=ones(N,N);
%          C3=ones(N,N);
%          C4=ones(N,N);
%         epsilon=2*epsilon;
%         remainder=rem(N,epsilon);
% 
%         for i=2:2*epsilon:N-remainder
%             for j=2:2*epsilon:N-remainder
%                 p=rand;
%                 if p<=1/2
%                    
%                
%                 C(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                  end
%             end
%         end
%         for i=epsilon+2:2*epsilon:N-remainder
%             for j=epsilon+2:2*epsilon:N-1-remainder
%                 p=rand;
%                 if p<=1/2
%                   
%                
%                 C(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                  end
%             end
%         end
%          for i=2:2*epsilon:N-remainder
%             for j=2+epsilon:2*epsilon:N-remainder
%                 p=rand;
%                 if p<=1/2
%                     
%                 
%                 C2(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                 end
%             end
%         end
%         for i=epsilon+2:2*epsilon:N-remainder
%             for j=2:2*epsilon:N-1-remainder
%                 p=rand;
%                 if p<=1/2
%                    
%                 
%                 C2(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                 end
%             end
%         end
%         
%         for i=2:2*epsilon:N-remainder
%             for j=2+epsilon:2*epsilon:N-remainder
%                 p=rand;
%                 if p<=1/2
%                     
%                 
%                 C3(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                 end
%             end
%         end
%         for i=epsilon+2:2*epsilon:N-remainder
%             for j=2:2*epsilon:N-1-remainder
%                 p=rand;
%                 if p<=1/2
%                    
%                 
%                 C3(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                 end
%             end
%         end
%         for i=2:2*epsilon:N-remainder
%             for j=2+epsilon:2*epsilon:N-remainder
%                 p=rand;
%                 if p<=1/2
%                    
%                 
%                 C4(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                 end
%             end
%         end
%         for i=epsilon+2:2*epsilon:N-remainder
%             for j=2:2*epsilon:N-1-remainder
%                 p=rand;
%                 if p<=1/2
%                     
%                
%                 C4(i:i+epsilon-2,j:j+epsilon-2)=c0*ones(epsilon-1,epsilon-1);
%                 end
%             end
%         end
%         
%         C1=hto2h(C,N);
%         C1=flipud(C1);
%         
%         C2=hto2h(C2,N);
%         C2=fliplr(C2(:,2:end));
%         C2=flipud(C2);
%         
%         C3=hto2h(C3,N);
%         C3=flipud(C3);
%         C3=flipud(fliplr(C3(1:end-1,2:end)));
%         
%         C4=hto2h(C4,N);
%         C4=C4(2:end,:);
%         
%         
%         
%         
%           C=[C2,C1;C3,C4];
%           

        for i=1:epsilon:N-epsilon
            for j=1:epsilon:N-epsilon
                p=rand;
                if p<=1/2
                    C(i:i+epsilon,j:j+epsilon)=2;
                else
                    C(i:i+epsilon,j:j+epsilon)=1;
                end
            end
        end
        
        
        
    case 3 %squares
         for i=1:N
            for j=1:N
                if rem(i-1,epsilon)>0 && rem(j-1,epsilon)>0
                    C(i,j)=2;
                end
            end
         end
         
%          for i=1:N
%             for j=1:N
%                 if i==1
%                     if rem(i-1,epsilon)==0
%                         C(i:i+4,j)=1;
%                     end
%                 
%                 elseif j==1
%                     if rem(j-1,epsilon)==0
%                         C(i,j:j+4)=1;
%                     end
%                 
%                 elseif i==N
%                     if rem(i-1,epsilon)==0
%                         C(i-4:i,j)=1;
%                     end
%                 
%                 elseif j==N
%                     if rem(j-1,epsilon)==0
%                         C(i,j-4:j)=1;
%                     end
%                 elseif rem(i-1,epsilon)==0
%                     C(i-2:i+2,j)=1;
%                 elseif rem(j-1,epsilon)==0
%                     C(i,j-2:j+2)=1;
%                 
%                 end
%             end
%          end
          
end


                
        
        
        


end

