function invp = getgrediant(node,pnearpoints)

pnode=node(pnearpoints,:);
pnode(:,1)=pnode(:,1)-pnode(1,1);
pnode(:,2)=pnode(:,2)-pnode(1,2);

P=zeros(6,6);
P(:,1)=1;
P(:,2)=pnode(:,1);
P(:,3)=pnode(:,2);
P(:,4)=pnode(:,1).^2;
P(:,5)=pnode(:,2).^2;
P(:,6)=pnode(:,1).*pnode(:,2);

invp=inv(P);

end