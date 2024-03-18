function K=AssmenbleK3D(K,invp,nid,nearpointsn)

    for i=1:size(nearpointsn,2)
        othernode=nearpointsn(1,i);
        K(nid,othernode)= K(nid,othernode)+2*invp(5,i)+2*invp(6,i)+2*invp(7,i);
    end
    
end