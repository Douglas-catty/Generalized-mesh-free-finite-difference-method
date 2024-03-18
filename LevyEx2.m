clear
clc

xline=linspace(-1,1,71);
%xline=linspace(-1,1,51);
points=[];
boundarypoints=[];
for i=1:size(xline,2)
  for j=1:size(xline,2)
    xx=xline(1,i);
    yy=xline(1,j);
    if xx>-1 & xx<1 &  yy>-1 & yy<1
    %if (xx^2+yy^2<1)
      points(:,end+1)=[xx;yy];
    else
      boundarypoints(:,end+1)=[xx;yy];  
    end
  end
end

%scatter(points(1,:),points(2,:));

h=xline(1,2)-xline(1,1);
alpha1=1;
cnalpha1=alpha1*gamma((1+alpha1)/2)/(2^(1-alpha1)*(pi^(1/2))*gamma(1-alpha1/2));
delta1=1;
alpha2=1;
cnalpha2=alpha2*gamma((1+alpha2)/2)/(2^(1-alpha2)*(pi^(1/2))*gamma(1-alpha2/2));
delta2=1;

KU=zeros(size(points,2),size(points,2));
FU=-1*ones(size(points,2),1);

for i=1:size(points,2)

  [i,size(points,2)]

  nowpoint=points(:,i);
  xx=nowpoint(1,1);
  yy=nowpoint(2,1);
  nowindexx=find(points(1,:)==xx); % points in the same column
  nowlinex=points(:,nowindexx);
  %scatter(nowlinex(1,:),nowlinex(2,:));hold on;
  a=min(nowlinex(2,:))-0.00001;
  b=max(nowlinex(2,:))+0.00001;
  nowindexy=find(points(2,:)==nowpoint(2,1)); % points in the same row
  nowliney=points(:,nowindexy);
  c=min(nowliney(1,:))-0.00001;
  d=max(nowliney(1,:))+0.00001;

  % points in the same column
  KU(i,i)=KU(i,i)-(delta2^alpha2)*(cnalpha2/alpha2)*(1/(yy-c)^alpha2+1/(d-yy)^alpha2);

  for j=1:size(nowlinex,2)
    if yy~=nowlinex(2,j)
      
      distancey=abs(yy-nowlinex(2,j));
      KU(i,nowindexx(1,j))=KU(i,nowindexx(1,j))+(delta2^alpha2)*cnalpha2*h/(distancey^(1+alpha2));
      KU(i,i)=KU(i,i)-(delta2^alpha2)*cnalpha2*h/(distancey^(1+alpha2));

    end
  end

  % points in the same row
  KU(i,i)=KU(i,i)+(delta1^alpha1)*cnalpha1*(-1)*(1/alpha1)*(1/((xx-a)^alpha1)+1/((b-xx)^alpha1));

  for j=1:size(nowliney,2)
    if xx~=nowliney(1,j)
      
      distancex=abs(xx-nowliney(1,j));
      KU(i,nowindexy(1,j))=KU(i,nowindexy(1,j))+(delta1^alpha1)*cnalpha1*h/(distancex^(1+alpha1));
      KU(i,i)=KU(i,i)-(delta1^alpha1)*cnalpha1*h/(distancex^(1+alpha1));

    end
  end

end

U=inv(KU)*FU;

scatter(points(1,:),points(2,:),[],U,'filled');
colormap('jet')
colorbar
axis equal
axis([-1.2,1.2,-1.2,1.2]);
title('Mean first exit time');

%Ualpha05=U;
%save('Ualpha05.mat','Ualpha05')

pause;

%%%%%%%
%Nboundary=201;
%theta=linspace(0,2*pi,Nboundary);
%theta=theta(1,1:Nboundary-1);
%costheta=cos(theta);
%sintheta=sin(theta);

%routside=1;

%boundarypoints=routside*[costheta;sintheta];

node=[points,boundarypoints];

%scatter(node(1,:),node(2,:))

num=size(node,2);
dis=zeros(num,num);%不同标号节点距离

for i=1:num
  for j=i+1:num
    
      dis(i,j)=norm(node(:,i)-node(:,j));
      dis(j,i)=dis(i,j);

  end
end

Taylor=6;
nearpoints=zeros(num,Taylor);
deltaB1=1;
deltaB2=1;
%6阶邻近点局部排序
for i=1:num
  for j=1:6
    [~,minnum]=min(dis(i,:));
    nearpoints(i,j)=minnum;
    dis(i,minnum)=max(dis(i,:));
  end
end

K2=zeros(num,num);% index2：Mean first exit time
F2=-1*ones(num,1);
for i=1:num
  xnode=node(1,i);
  ynode=node(2,i);
  %r=norm([xnode,ynode]);

  %if abs(r-1)<1e-5 %边界
   if (xnode==-1) | (xnode==1) | (ynode==-1) | (ynode==1)
        K2(i,i)=1;
        F2(i)=0;
  
  else
     %内部点
     invp=getgrediant(node,nearpoints(i,:));
     K2=AssmenbleK(K2,invp,i,nearpoints(i,:),deltaB1,deltaB2);
  end



end

U2=inv(K2)*F2;

scatter(node(1,:),node(2,:),[],U2,"filled")
colormap('jet')
colorbar
axis equal
axis([-1.5,1.5,-1.5,1.5]);
title('Mean first exit time');
pause;

K2base=K2;
K2base(1:size(points,2),1:size(points,2))=K2base(1:size(points,2),1:size(points,2))+KU;
U2base=inv(K2base)*F2;
scatter(node(1,:),node(2,:),[],U2base,"filled")
colormap('jet')
colorbar
axis equal
axis([-1.5,1.5,-1.5,1.5]);
title('Mean first exit time');
