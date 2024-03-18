clear
clc

load node.mat
r=0.5;
R=1;

num=size(node,2);
dis=zeros(num,num);%不同标号节点距离

for i=1:num
    [i,num]
  for j=i+1:num
    
      dis(i,j)=norm(node(:,i)-node(:,j));
      dis(j,i)=dis(i,j);

  end
end

%Taylor=10;
Taylor=19;
nearpoints=zeros(num,Taylor);
%10阶邻近点局部排序
for i=1:num
    [i,num]
  for j=1:Taylor
    [~,minnum]=min(dis(i,:));
    nearpoints(i,j)=minnum;
    dis(i,minnum)=max(dis(i,:));
  end
end

%重构矩阵K
K1=zeros(num,num);% index1：Exit probability
F1=zeros(num,1);
K2=zeros(num,num);% index2：Mean first exit time
F2=-1*ones(num,1);

for i=1:num

  [i,num]

  xnode=node(1,i);
  ynode=node(2,i);
  znode=node(3,i);

  normnode=norm([xnode,ynode,znode]);

  %if abs(normnode-r)<1e-4 %内边界
  if abs(normnode-r)<0.02 %内边界
        K1(i,i)=1;
        F1(i)=1;
        %K2(i,i)=1;
        F2(i)=0;
  %elseif abs(normnode-R)<1e-4 %外边界
  elseif abs(normnode-R)<0.02 %外边界
        K1(i,i)=1;
        F1(i)=0;
        %K2(i,i)=1;
        F2(i)=0;

  else
     %内部点
     invp=getgrediant3D(Taylor,node,nearpoints(i,:));
     K1=AssmenbleK3D(K1,invp,i,nearpoints(i,:));
     K2=K1;
  end



end

U1=inv(K1)*F1;
U2=inv(K2)*F2;

figure
set(0,'defaultfigurecolor','w')
scatter3(node(1,:),node(2,:),node(3,:),6*(abs(U1)+0.1),U1)
colormap('jet')
colorbar
axis equal
axis([-1.2,1.2,-1.2,1.2]);
title('Escape probability');
xlabel('X');
ylabel('Y');
zlabel('Z');

