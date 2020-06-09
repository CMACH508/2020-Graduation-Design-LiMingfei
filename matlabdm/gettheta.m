function theta=gettheta(vec, row, col, beta, hv)
    rr = beta*log(row);
    cc = beta*log(col);
    ksai = 0.6;gamma=0.8;m = 0;
    while m<1000
        v = vec + ksai^m * (hv - vec);
        L= 0;
        for i =1:n
            for j =1:n
                for k=1:n-1
                    L=L+d[i][j] * v[i * n + k] * v[j * n + k + 1]
                end
            L =L+ d[i][j] * v[i * n + k] * v[j * n + n - 1]- 0.5 * rho * v[i * n + j]^2;
            end
        end
        tmp = 0;
        for i =1:n
            for j =1:n
                tmp =tmp+v(i+(j-1)*n);
            end
            tmp=tmp-1;
            L=L+rr(i)*tmp+cc(i)*tmp;
        end
        L2 = e(vec);
        for i=1:n
            for j =1:n
                tmp=tmp+vec(i+(j-1)*n);
            end
            tmp=tmp-1;
            L2=L2+rr(i)*tmp+cc(i)*tmp;
        end
        tmp1 = gamma*(hv - vec);
        tmp2 = gradient(vec, rr, cc)';
        L2 =L2+ ksai^m*tmp1*tmp2;
        if L <= L2
            break
        end
        m=m+1;
    end
    theta=ksai^m;
end