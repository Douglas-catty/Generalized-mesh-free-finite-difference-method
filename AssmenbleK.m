function K=AssmenbleK(K,invp,nid,pset)

    for i=1:length(pset)
        othernode=pset(i);
        K(nid,othernode)= K(nid,othernode)+2*invp(4,i)+2*invp(5,i);
    end
    
end