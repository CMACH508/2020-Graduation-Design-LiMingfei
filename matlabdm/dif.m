function d=dif(dist,vec,n,i,rho)
        e=0;
        if mod(i,n)~=0 
            r=mod(i,n);
        else r=n;
        end       
        c=ceil(i/n);
        d=0;
        if c==1
            for k=1:n
                e= e+dist(k,r)*vec(k+n^2-n)+dist(r,k)*vec(k+n);
                if k+n^2-n==i 
                    d=d+dist(k,r);
                end
                if k+n==i
                    d=d+dist(r,k);
                end
            end
        end
        if c==n
           for k=1:n
                e= e+dist(k,r)*vec(k+n^2-2*n)+dist(r,k)*vec(k);
                if k+n^2-2*n==i 
                    d=d+dist(k,r);
                end
                if k==i
                    d=d+dist(r,k);
                end
           end    
        end
        if c>1 && c<n
           for k=1:n
                e= e+dist(k,r)*vec(k+(c-2)*n)+dist(r,k)*vec(k+c*n);
                if k+(c-2)*n==i 
                    d=d+dist(k,r);
                end
                if k+c*n==i
                    d=d+dist(r,k);
                end
           end
        end
        d=d+e-rho*vec(i);
end