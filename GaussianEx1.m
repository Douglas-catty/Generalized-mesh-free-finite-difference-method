%%% 生成边界结点（边界节点越密越好？）
%fd=@(p) ddiff(p(:,1).^2+p(:,2).^2-9,p(:,1).^2+p(:,2).^2-1);
%[p,t]=distmesh2d(fd,@huniform,0.21,[-3,-3;3,3],[]);
%pause

%%% 域内节点均匀分布
%%% boundary points : as dense as possible
Nboundary=201;
theta=linspace(0,2*pi,Nboundary);
theta=theta(1,1:Nboundary-1);
costheta=cos(theta);
sintheta=sin(theta);

rinside=1;
routside=3;

boundaryinside=rinside*[costheta;sintheta]';
boundaryoutside=routside*[costheta;sintheta]';
boundarypoints=[boundaryinside;boundaryoutside];

%%% inner points are uniformiy distributed
%%% uniformly distributed in interval [a,b]:a+(b-a)*rand(n,m)
%%% Halton序列生成均匀分布
ax=-3.5;
bx=3.5;
ay=-3.5;
by=3.5;
Nx=50+1;% +1保证（0，0）位于格点上
Ny=50+1;
Nxyindex=[];
m=1;

for i=1:Nx
    for j=1:Ny
      Nxyindex(:,end+1)=[i;j];
      m=m+1;
    end
end

uniformpointsx=ax+(bx-ax)*(Nxyindex(1,:)./Nx);
uniformpointsy=ay+(by-ay)*(Nxyindex(2,:)./Ny);
uniformpoints=[uniformpointsx',uniformpointsy'];

Runiformpoints=sqrt(uniformpoints(:,1).^2+uniformpoints(:,2).^2);

innerpoints=uniformpoints(find(Runiformpoints>1 & Runiformpoints<3),:);

node=[innerpoints;boundarypoints];

scatter(node(:,1),node(:,2));





num=size(node,1);
dis=zeros(num,num);%不同标号节点距离

for i=1:num
  for j=i+1:num
    
      dis(i,j)=norm(node(i,:)-node(j,:));
      dis(j,i)=dis(i,j);

  end
end

Taylor=6;
nearpoints=zeros(num,Taylor);
%6阶邻近点局部排序
for i=1:num
  for j=1:6
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
  xnode=node(i,1);
  ynode=node(i,2);
  r=norm([xnode,ynode]);

  if abs(r-1)<1e-5 %内边界
        K1(i,i)=1;
        F1(i)=1;
        K2(i,i)=1;
        F2(i)=0;
  elseif abs(r-3)<1e-6 %外边界
        K1(i,i)=1;
        F1(i)=0;
        K2(i,i)=1;
        F2(i)=0;

  else
     %内部点
     invp=getgrediant(node,nearpoints(i,:));
     K1=AssmenbleK(K1,invp,i,nearpoints(i,:));
     K2=AssmenbleK(K1,invp,i,nearpoints(i,:));
     F1(i)=0;
  end



end

U1=inv(K1)*F1;
U2=inv(K2)*F2;

subplot(1,2,1);
set(0,'defaultfigurecolor','w')
scatter(node(:,1),node(:,2),[],U1,"filled")
colormap('jet')
colorbar
axis equal
axis([-3.5,3.5,-3.5,3.5]);
title('Escape probability');

subplot(1,2,2);
scatter(node(:,1),node(:,2),[],U2,"filled")
colormap('jet')
colorbar
axis equal
axis([-3.5,3.5,-3.5,3.5]);
title('Mean first exit time');