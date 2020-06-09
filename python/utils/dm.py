# DM
import numpy as np
import math
import sys
import datetime
epsilon = 0.01
beta_0 = 200
mu_k = 0.95
rho = 80

def de_dvij(d,i,j,v,rho,n):
    temp = 0
    if j == 0:
        for k in range(n):
            temp += d[k, i] * v[k, n - 1] + d[i, k] * v[k, j + 1]
    elif j == n - 1:
        for k in range(n):
            temp += d[k, i] * v[k, j - 1] + d[i, k] * v[k, 0]
    else:
        for k in range(n):
            temp += d[k, i] * v[k, j - 1] + d[i, k] * v[k, j + 1]
    temp-=rho*v[i,j]
    return temp

def alpha_ij_v(d,i, j, v, beta,rho,n):
    temp=math.exp(de_dvij(d,i,j,v,rho,n)/beta)
    return temp
def xi_rc(d,i, r, c, v, beta,rho,n):
    temp = 0
    for j in range(n):
        alpha = alpha_ij_v(d,i, j, v, beta,rho,n)
        temp += 1 / (1 + r[i] * c[j] * alpha)
    temp = temp - 1
    return temp * r[i]


def yj_rc(d,j, r, c, v, beta,rho,n):
    temp = 0
    for i in range(n):
        alpha = alpha_ij_v(d,i, j, v, beta,rho)
        temp += 1 / (1 + r[i] * c[j] * alpha)
    temp = temp - 1
    return temp * c[j]
def f(d,r, c, v, beta,rho,n):
    result = 0
    for i in range(n):
        part1 = 0
        for j in range(n):
            alpha = alpha_ij_v(d,i, j, v, beta,rho,n)
            part1 += 1 / (1 + r[i] * c[j] * alpha)
        part1 -= 1
        part1 = part1 ** 2
        result += part1
    for j in range(n):
        part2 = 0
        for i in range(n):
            alpha = alpha_ij_v(d,i, j, v, beta,rho,n)
            part2 += 1 / (1 + r[i] * c[j] * alpha)
        part2 -= 1
        part2 = part2 ** 2
        result += part2
    result = result / 2
    return result
# fix v
def iterate(d,v, r, c, beta,rho,n):
    x = np.zeros(n)
    y = np.zeros(n)
    while np.sqrt(f(d,r, c, v, beta,rho,n)) >= 0.001:
        # initialize x(r,c) y(r,c)
        for i in range(n):
            x[i] = xi_rc(d,i, r, c, v, beta,rho,n)
        for j in range(n):
            y[j] = yj_rc(d,j, r, c, v, beta,rho,n)
        r = r + mu_k * x
        c = c + mu_k * y
    return r, c

# fix r,c
def calculate_v(d,v,r, c,beta,rho,n):
    for i in range(n):
        for j in range(n):
            alpha = alpha_ij_v(d,i, j, v, beta,rho,n)
            v[i , j] = 1 / (1 + r[i] * c[j] * alpha)
    return v
def e(d,v,beta,rho,n):
    result = 0
    for i in range(n):
        for j in range(n):
            for k in range(n - 1):
                result += d[i,j] * v[i, k] * v[j,k + 1]
            result += d[i,j] * v[i ,n-1] * v[j,0]
            result -= 0.5 * rho * v[i, j] ** 2
    temp=0
    for i in range(n):
        for j in range(n):
            if v[i,j] != 0 and v[i , j] != 1:
                temp+=v[i,j]*np.log(v[i,j])+(1-v[i,j])*np.log(1-v[i,j])
    result+=beta*temp
    return result
def gradient(d,v, lamda_rk, lamda_ck,beta,n):
    dL = np.zeros_like(v)
    for i in range(n):
        for j in range(n):
            if j == 0:
                for k in range(n):
                    dL[i, j] = dL[i , j] + d[k,i] * v[k , n - 1] + d[i,k] * v[k , j + 1]
            elif j == n - 1:
                for k in range(n):
                    dL[i ,j] = dL[i,j] + d[k,i] * v[k ,j - 1] + d[i,k] * v[k,0]
            else:
                for k in range(n):
                    dL[i,j] = dL[i,j] + d[k,i] * v[k , j - 1] + d[i,k] * v[k,j + 1]
            dL[i,j] -= rho * v[i, j]
            dL[i ,j] = dL[i, j] + lamda_rk[i] + lamda_ck[j] + beta * np.log(v[i, j] / (1 - v[i, j]))
    return dL



def theta_k(d,v, r, c, beta, h,rho,n):
    lamda_rk = beta * np.log(r)
    lamda_ck = beta * np.log(c)
    xi = 0.6
    gamma = 0.8
    mk = 0
    while True:
        v1 = v + math.pow(xi, mk) * (h - v)
        L1 = e(d,v1,beta,rho,n)
        temp = 0
        for i in range(n):
            for j in range(n):
                temp += v1[i,j]
            temp -= 1
            L1+=lamda_rk[i] * temp
        temp=0
        for j in range(n):
            for i in range(n):
                temp+=v1[i,j]
            temp-=1
            L1+=lamda_ck[j]*temp

        L2 = e(d,v,beta,rho,n)
        temp=0
        for i in range(n):
            for j in range(n):
                temp += v[i, j]
            temp -= 1
            L2 = L2 + lamda_rk[i] * temp
        temp=0
        for j in range(n):
            for i in range(n):
                temp+=v[i,j]
            temp-=1
            L2+=lamda_ck[j]*temp


        temp1 = gamma * (h - v).reshape(1, n * n)
        temp2 = gradient(d,v, lamda_rk, lamda_ck,beta,n).reshape(n * n, 1)

        L2 += math.pow(xi, mk) * np.dot(temp1, temp2)[0, 0]
        if L1 <= L2:
            result = math.pow(xi, mk)
            break
        else:
            mk = mk + 1
        # print("mk", mk)
        if mk>5:
            mk=5
            result = 0.1
            break


    return result
def translate_v(v):
    newv=np.zeros_like(v)
    for i in range(n):
        for j in range(n):
            if v[i,j]>=0.9:
                newv[i,j]=1
    return newv
def belong_P(newv):
    for i in range(n):
        for j in range(n):
            if newv[i,j]>1 or newv[i,j]<0:
                return False
    temp=newv[:].sum(axis=0)
    for i in range(n):
        if temp[i]!=1:
            return False
    temp=newv[:].sum(axis=1)
    for i in range(n):
        if temp[i]!=1:
            return False
    return True




def dmethod(y_preds, x_edges_values, batch_size, num_nodes)
    start = datetime.datetime.now()
    # d = np.loadtxt('city_30.txt')
    d=x_edges_values
    n=num_nodes
    #v_ = np.random.rand(n,n)
    v_=y_preds
    print(v_)
    r_0 = np.ones(n)
    c_0 = np.ones(n)
    #step 0
    v = v_
    r = r_0
    c = c_0
    v_q = v
    beta = beta_0
    r, c = iterate(d,v, r, c, beta,rho,n)
    new_v = calculate_v(d,v, r, c, beta,rho,n)
    #step 1
    while True:
        v=new_v
        #print("here")
        r, c = iterate(d,v, r, c, beta,rho,n)
        h=np.zeros_like(v)
        for i in range(n):
            for j in range(n):
                alpha=alpha_ij_v(d,i,j,v,beta,rho,n)
                h[i,j]=1/(1+r[i]*c[j]*alpha)
        print("beta",beta)
        if np.linalg.norm(h - v, ord=2) < epsilon:
            if beta < 1:
                print("process1 terminate")
                break
            else:
                v_q = v
                new_v=v
                beta = 0.95 * beta
        else:
            theta = theta_k(d,v, r, c, beta, h,rho,n)
            print("theta",theta)
            new_v = v + theta * (h - v)
    # step 0
    beta = 1
    v = v_q
    new_v = v
    #step 1
    v_t = translate_v(v[:])
    if belong_P(v_t):
        #print("process2 terminate without enter")
        #print("d", d)
        dis = 0
        for i in range(n):
            for j in range(n):
                if v_t[i, j] == 1:
                    dis += d[i][j]
        #print("dis", dis)
        #print("v", v_t)
        end = datetime.datetime.now()
        #print(end - start)
    else:
        rho+=2
    while True:
        #step2
        v=new_v
        r, c = iterate(d,v, r, c, beta,rho,n)
        #step3
        h = np.zeros_like(v)
        for i in range(n):
            for j in range(n):
                alpha = alpha_ij_v(d,i, j, v, beta,rho,n)
                h[i, j] = 1 / (1 + r[i] * c[j] * alpha)
        if np.linalg.norm(h - v, ord=2) < epsilon:
            # step1
            #print("Step1")
            #print("rho",rho)
            v_t=translate_v(v[:])
            if belong_P(v_t):
                #print("process2 terminate")
                #print("d", d)
                dis = 0
                for i in range(n):
                    for j in range(n):
                        if v_t[i, j] == 1:
                            dis += d[i,j]
                #print("dis", dis)
                #print("v", v_t)
                end = datetime.datetime.now()
            else:
                rho += 2
        else:
            theta = theta_k(d,v, r, c, beta, h,rho,n)
            new_v = v + theta * (h - v)
            
    return v