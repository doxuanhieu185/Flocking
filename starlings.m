
function [X1 X2 V1 V2 ]=starlings()
it=500;
n=30;
connectivityradius=1;
X1=zeros(it,n);
X2=zeros(it,n);
X1(1,:)=rand(1,n);
X2(1,:)=rand(1,n);

V1=zeros(it,n);
E1=zeros(it,n);
V2=zeros(it,n);
V21=zeros(it,n);
E21=zeros(it,n);
V22=zeros(it,n);
E22=zeros(it,n);


V1(1,:)=0.1*rand(1,n);
V2(1,:)=(2*pi)*rand(1,n);
Neighbors=repmat(struct('val',[]), n, 1);

epsilon=0.005; 
epsilons1=epsilon*ones(n,1);
counter1=ones(n,1);
alpha=0.1;
Silent1=false(n,1);

epsilons21=epsilon*ones(n,1);
counter21=ones(n,1);
Silent21=false(n,1);

epsilons22=epsilon*ones(n,1);
counter22=ones(n,1);
Silent22=false(n,1);

pass=0;
for kk=2:it
    [A I]=MakeGraph(X1(1,:),X2(1,:),connectivityradius);
    W=0.2*eye(n)+0.8*WldCreate(I);
    for j=1:n
     Neighbors(j).val=find(W(j,:)); 
    end
    
    for jj=1:3
        if jj==1
            x=V1(kk-1,:)';
            e=E1(kk-1,:)';
            epsilons=epsilons1;
            counter=counter1;
            Silent=Silent1;
        else
            if jj==2
                V21(kk-1,:)=sin(V2(kk-1,:));
                x=V21(kk-1,:)';
                e=E21(kk-1,:)';
                epsilons=epsilons21;
                counter=counter21;
                Silent=Silent21;
            else
                V22(kk-1,:)=cos(V2(kk-1,:));
                x=V22(kk-1,:)';
                e=E22(kk-1,:)';
                epsilons=epsilons22;
                counter=counter22;
                Silent=Silent22;
            end
        end
    y=W*x;   
    d=y+e-x;
    wait=[];
    waits=false(n,1);
    
    MaxMin=eye(n);
    transmit=[];
    xoriginal=x;
    
    for i=1:n        
        
        if abs(d(i))<epsilons(i)
            
            if ~Silent(i)
                   Silent(i)=true;
            end
            
            Ind=find(W(i,:));
            v=xoriginal(Ind);
            z1=0;
            z2=0;
            
            for j=1:length(Ind)
              if v(j)-xoriginal(i)>0
                z1=z1+W(i,Ind(j))*v(j);
                z2=z2+W(i,Ind(j))*xoriginal(i);
              else
                z1=z1+W(i,Ind(j))*xoriginal(i);
                z2=z2+W(i,Ind(j))*v(j);
              end
            end
            
            if  abs(y(i)+e(i)-z1)<epsilons(i) && abs(y(i)+e(i)-z2)<epsilons(i)
                wait(length(wait)+1)=i; 
                waits(i)=true;
               % e(i)=d(i);
            else
                
                  pass=1;
                  
            end
          
        end
        
        if abs(d(i))>=epsilons(i) || pass
            pass=0;
            if Silent(i)
                Silent(i)=false;
                
               if rand<0.5 
                 counter(i)=counter(i)+1;
                 epsilons(i)=epsilons(i)+epsilon/(counter(i)^2);
               end
            end
            
            CC=xoriginal(Neighbors(i).val);
            maxn=find(CC==max(CC));
            minn=find(CC==min(CC));
            MaxMin(Neighbors(i).val(maxn(1)),i)=1;
            MaxMin(Neighbors(i).val(minn(1)),i)=1;
            
            transmit(length(transmit)+1)=i;
            %M(k)=M(k)+1;
            
                c=alpha*(d(i)-e(i))/(1-W(i,i));
            if abs(e(i))>abs(c)
                
                x(i)=y(i)+sign(c*e(i))*c;
                e(i)=e(i)-sign(c*e(i))*c;
            else
                x(i)=y(i)+e(i);
                e(i)=0;
            end
        end
               
    end
 
  
   
  for j=1:length(wait)
      i=wait(j);
      Ind=Neighbors(i).val;
      %VV=waits(Ind);
      VV=MaxMin(i,Ind);
      if sum(VV)==1
          e(i)=d(i);
      else
          if Silent(i)
                Silent(i)=false;
                
               if rand<0.5
                   counter(i)=counter(i)+1;
                   epsilons(i)=epsilons(i)+epsilon/(counter(i)^2);
               end
          end
          
         % M(k)=M(k)+1;
          
            %VDiff=setdiff(Ind,wait);
            v=xoriginal(Ind);
            z=0;
            for l=1:length(Ind)
              %if ismember(Ind(l),VDiff)
              if VV(l)
                z=z+W(i,Ind(l))*v(l);
              else
                z=z+W(i,Ind(l))*xoriginal(i);
              end
            end
            
            e(i)=y(i)-z+e(i);
            
            x(i)=z;     
      end
  end
  
    if jj==1
            V1(kk,:)=x';
            E1(kk,:)=e';
            epsilons1=epsilons;
            counter1=counter;
            Silent1=Silent;
    else
        if jj==2
            V21(kk,:)=x';
            E21(kk,:)=e';
            epsilons21=epsilons;
            counter21=counter;
            Silent21=Silent;
        else
            V22(kk,:)=x';
            E22(kk,:)=e';
            epsilons22=epsilons;
            counter22=counter;
            Silent22=Silent;
        end
    end
    
    end
    V2(kk,:)=mod(atan2(V21(kk,:),V22(kk,:)),2*pi);
    X1(kk,:)=X1(kk-1,:)+V1(kk,:).*cos(V2(kk,:));
    X2(kk,:)=X2(kk-1,:)+V1(kk,:).*sin(V2(kk,:));
  % X1(kk,:)=X1(kk-1,:);
  if (any(X1(kk,:)<0))
    OUT=find(X1(kk,:)<0);
    V2(kk,OUT)=mod(pi-V2(kk,OUT),2*pi);
    X1(kk,OUT)=-X1(kk,OUT);
  end
  
  if (any(X1(kk,:)>1))
    OUT=find(X1(kk,:)>1);
    V2(kk,OUT)=mod(pi-V2(kk,OUT),2*pi);
    X1(kk,OUT)=2-X1(kk,OUT);
  end
    
  
    
    if (any(X2(kk,:)<0))
    %X2(kk,:)=X2(kk-1,:);
    OUT=find(X2(kk,:)<0);
    V2(kk,OUT)=mod(-V2(kk,OUT),2*pi);
    X2(kk,OUT)=-X2(kk,OUT);
  end
  
  if (any(X2(kk,:)>1))
    OUT=find(X2(kk,:)>1);
    V2(kk,OUT)=mod(-V2(kk,OUT),2*pi);
    X2(kk,OUT)=2-X2(kk,OUT);
  end
    
    
    
end





plotstarlings(X1,X2,it);

end

function [A I]=MakeGraph(x,y,c)
    
n=length(x);
z=zeros(n,2);
z(:,1)=x;
z(:,2)=y;
A=zeros(n,n);
I=zeros(n,1);
s=0;
r=sqrt(c*log(n)/n);
for i=1:n
    for j=(i+1):n
        if (sqrt((z(i,1)-z(j,1))^2+ (z(i,2)-z(j,2))^2)<=r)        
            A(i,j)=1;
            A(j,i)=1;
            s=s+1;
            I(i,s)=1;
            I(j,s)=-1;
        end
    end
end
end

function plotstarlings(X1,X2,K)
figure;
 for it=1:K-1
hLine=plot(X1(it:1+it,:),X2(it:it+1,:),'LineWidth',2);
axis([0 1 0 1])
pause(0.2);
if ~ishandle(hLine)
    break;
end
end

end

function Wld=WldCreate(I)
[l h]=size(I);
L=zeros(l,l);
 for i=1:h
    j=find(I(:,i)>0);
    k=find(I(:,i)<0);
    L(j,k)=1;
 end
L=L+L';
D=zeros(l,1);
Wld=zeros(l,l);
D=sum(L);

 for i=1:l
    for j=1:l
        if (L(i,j)==1) Wld(i,j)=1/(1+max(D(i),D(j)));
        end
    end
 end

 for i=1:l
     Wld(i,i)=1-sum(Wld(:,i));
 end
end