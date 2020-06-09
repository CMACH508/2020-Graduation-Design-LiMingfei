function v = getv(dist,vec,rol,col,beta,n,rho)   
    for i=1:n^2
        if mod(i,n)~=0 
            rr=rol(mod(i,n));
        else rr=rol(n);
        end
        cc=col(ceil(i/n));        
        d=dif(dist,vec,n,i,rho);    
        v(i)=1/(exp(rr/beta)*exp(cc/beta)*exp(d/beta));
    end
end