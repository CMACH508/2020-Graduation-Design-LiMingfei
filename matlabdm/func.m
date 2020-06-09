function f=func(dist,row,col,vec,beta,n,rho)
    f=0;
    for i =1:n
        tmp = 0;
        for j =1:n
            cnt=j*n+i;
            alpha = exp(dif(dist,vec,n,cnt,rho)/beta);
            tmp =tmp+ 1/(1 + row(i)* col(j) * alpha);
        end
        tmp=(tmp-1)^2;
        f = f+tmp;
    end
    for j =1:n
        tmp = 0;
        for i =1:n
            cnt=j*n+i;
            alpha = exp(dif(dist,vec,n,cnt,rho)/beta);
            tmp =tmp+ 1/(1 + row(i)* col(j) * alpha);
        end
        tmp=(tmp-1)^2;   
        f=f+tmp;
    end
    f=sqrt(f/2);
end