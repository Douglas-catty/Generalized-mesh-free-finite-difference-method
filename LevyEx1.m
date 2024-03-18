clear;
clc;

N=2001;
meshx=linspace(-2,2,N);
h=meshx(1,2)-meshx(1,1);
alpha=0.5;
cnalpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*(pi^(1/2))*gamma(1-alpha/2));

K=zeros(N,N);
F=zeros(N,1);% F for mean exit time
F2=zeros(N,1);% F2 for exit probability distribution


for i=1:501

    K(i,i)=1;

end

for i=1501:2001

    K(i,i)=1;
    F2(i,1)=1;

end

for i=502:1500

      F(i,1)=-1;
      K(i,i)=-(cnalpha/alpha)*(1/(1+meshx(1,i))^alpha+1/(1-meshx(1,i))^alpha);

  for k=501-i:1:1501-i

      if k~=0
        K(i,i+k)=K(i,i+k)+h*cnalpha/(abs(meshx(1,i+k)-meshx(1,i))^(1+alpha));
        K(i,i)=K(i,i)-1*h*cnalpha/(abs(meshx(1,i+k)-meshx(1,i))^(1+alpha));
      end

  end
        
        F2(i,1)=-1*cnalpha*(1/alpha)*(1/(1-meshx(1,i))^alpha);

end

U=inv(K)*F;
P=inv(K)*F2;

%%% 理论解

ureal=@(x) (sqrt(pi)*(1^2-x.^2).^(alpha/2))./((2^alpha)*gamma(1+alpha/2)*gamma(1/2+alpha/2));
Ureal=ureal(meshx);

ExitPIntegral=@(x) (1-x.^2).^(alpha/2-1);
interval=meshx(1,501:1501);
Pcoefficient=((2*1)^(1-alpha))*gamma(alpha)/(gamma(alpha/2)^2);
Preal=0*interval;

for i=1:size(interval,2)

Preal(1,i)=Pcoefficient*quad(ExitPIntegral,-1,interval(1,i),0.0000001);

end

% symmetry of exit probability：P（x）+P（-x）=1
for i=1:500

    Preal(1,501-i)=1-Preal(1,501+i);

end

%%%

subplot(1,2,1);
set(0,'defaultfigurecolor','w')

plot(meshx,U,'g-','linewidth',2);
hold on;
plot(meshx,Ureal,'r--','linewidth',2);
axis([-1,1,0,(1.2)*max(U)]);

legend('Numerical solution','Analytical solution');
title('Mean first exit time');

subplot(1,2,2);
plot(meshx,P,'g-','linewidth',2);
hold on;
plot(interval,Preal,'r--','linewidth',2);
axis([-1,1,0,1.2]);
legend('Numerical solution','Analytical solution');
title('Escape probability');

%save('alpha05','meshx','U','Ureal','P','interval','Preal')