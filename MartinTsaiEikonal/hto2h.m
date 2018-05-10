function [ U2h] = hto2h( Uh,Nh )
N=((Nh-1)/2)+1;
U2h=zeros(N,N);
for i=1:N
    for j=1:N
       U2h(i,j)=Uh(2*i-1,2*j-1); 
    end
end
end
