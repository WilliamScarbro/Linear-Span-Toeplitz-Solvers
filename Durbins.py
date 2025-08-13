from Util import *
from StabilityUtil import *
from functools import reduce

def durbins_prog(p,k,y,b,ops):
    alpha = - (p[k] + dot(p[1:k],rev(y))) / b # (p[0] + dot(p[1:k],y))
    b=(1-alpha**2)*b
    # print("alpha",alpha)
    y_new = [y[i]+alpha*y[k-2-i] for i in range(k-1)]+[alpha]
    ops=(ops[0]+2*k,ops[1]+2*k)
    return y_new,b,ops

def durbins(p,N):
    # print(p)
    # print(N)
    y=[-p[1]/p[0]]
    b=p[0] + p[1]*y[0]

    ops=(0,0)
    for k in range(2,N+1):
        y,b,ops=durbins_prog(p,k,y,b,ops)
        
    return y,ops

def durbins_linear_span_prog(N,p,k,Y,Y_hat,y,ops,for_levinsons=False):
    alpha = - (p[k] + Y_hat[0]) / (p[0] + Y[0])
    # true_alpha = - (p[k] + dot(p[1:k],rev(y))) / (p[0] + dot(p[1:k],y))
    
    Y_index_max = N-k
    Y_hat_new=[Y_hat[q+1]+alpha*(p[1+q]+Y[q+1]) for q in range(0,Y_index_max)]
    # print("check",Y[0],alpha,p[k],Y_hat[0])
    Y_new=[Y[q]+alpha*(p[k+q]+Y_hat[q]) for q in range(0,Y_index_max)]
    y_new = [y[i]+alpha*y[k-2-i] for i in range(k-1)]+[alpha]


    # print("alpha diff",alpha,true_alpha-alpha)
    # print("Y_hat diff",[Y_hat_new[q]-dot(p[1+q:k+q+1],rev(y_new)) for q in range(N-k-1)])
    # print("Y true",[dot(p[1+q:k+q+1],y_new) for q in range(N-k-1)])
    # print("Y_new",Y_new)
    ops=(ops[0]+4*(N-k)+k,ops[1]+2*(N-k)+k)
    
    return Y_new,Y_hat_new,y_new,ops
    
def durbins_linear_span(p,N):
    y=[-p[1]/p[0]]
    Y_hat=[p[1+q]*y[0] for q in range(N)]
    Y=Y_hat.copy()
    # print(p)
    # print(y)
    # print(Y_hat)
    # print(Y)

    ops=(0,0)
    for k in range(2,N+1):
        Y,Y_hat,y,ops = durbins_linear_span_prog(N,p,k,Y,Y_hat,y,ops)
        # alpha = - (p[k] + r[0]) / (p[0] + d[0])
        # r_new=[r[q+1]+alpha*(p[1+q]-d[q+1]) for q in range(0,N-k+1)]
        # d_new=[d[q]+alpha*(p[k+q]-r[q]) for q in range(0,N-k+1)]
        # r=r_new
        # d=d_new
        # y_new = [y[i]+alpha*y[k-2-i] for i in range(k-1)]+[alpha]
        # y=y_new
        
    return y,ops

# 
def durbins_split_prog(p,k,u,U,ops):
    # print("prog",p,k,u,U)
    # print(p[0],p[k],p[k-1])
    gamma = (p[0] + U[1] + p[k]) / (p[0] + U[0] + p[k-1])
    # gamma check
    # print("gamma*",gamma)
    # print("gamma(0)",u[1][0] + 1 + (p[1] + p[k]) )
    # for i in range(1,k-1):
    #     print(f"gamma({i})",(u[1][i] + u[1][i-1] -(p[1+i]+p[k-i])/p[i])/u[0][i-1])
    # print("alpha",alpha)
    # print("gamma(alpha)",(1-alpha[1])*(1+alpha[0]))
    
    u_new_len = (k-1)//2
    # print(range(1,(k-1)//2+1))
    u_new = [u[1][0]-gamma+1]+[u[1][i] + u[1][i-1] - gamma*u[0][i-1] for i in range(1,(k-1)//2)]
    # for odd, looking back at last element misses end of half-list (handle as separate case)
    if k%2==1:
        u_new.append(u[1][-1] + u[1][-1] - gamma*u[0][-1])
    else:
        u_new.append(u[1][-1] + u[1][-2] - gamma*u[0][-1])
        
    # print("u_new",u_new)
    if k%2==0:
        # print(p[1:k//2+1],p[k//2+1:k+1])
        # print(u_new,rev(u_new))
        U_new = dot(p[1:k//2+1],u_new) + dot(p[k//2+1:k+1],rev(u_new))
    else:
        # print(p[1:k//2+2],p[k//2+2:k+1])
        # print(u_new,rev(u_new[:-1]))
        U_new = dot(p[1:k//2+2],u_new) + dot(p[k//2+2:k+1],rev(u_new[:-1]))

    # print("U_new",U_new)
    ops=(ops[0]+(k-1)+k,ops[1]+(k-1)//2+k//2)
    return (u[1],u_new),(U[1],U_new),ops


def durbins_split(p,N):
    # print("p",p)
    # dubrin init
    y0=[-p[1]/p[0]]
    b0=p[0] + p[1]*y0[0]

    y1,b1,_ = durbins_prog(p,2,y0,b0,(0,0))
    y=(y0,y1)

    y_check=y1
    beta_check=b1
    
    # print(y)
    # we store all of u to make U initialization easier. Later, we only store first half of u (u is symmetric)
    u=tuple([y[i][j]+y[i][i-j] for j in range(i+1)] for i in range(2))
    # print("u",u)
    U=tuple(dot(u[i],p[1:i+2]) for i in range(2))
    # print("U",U)



    #debug
    # alpha = (get_alpha(2,u),get_alpha(1,(
    
    # fix u so only half is stored
    u=(u[0],[u[1][0]])
    ops=(0,0)
    for k in range(3,N+1):
        # print("k",k)
        u,U,ops = durbins_split_prog(p,k,u,U,ops)

        
    # recover y
    sum_from_half = lambda n,l : 2*sum(l) if n%2==0 else 2*sum(l)-l[-1]
    sum_u_n = sum_from_half(N,u[1])
    sum_u_nmo = sum_from_half(N-1,u[0])
    alpha = (sum_u_n - sum_u_nmo) / (sum_u_nmo + 2)

    u=(fill_from_half(N-1,u[0]),fill_from_half(N,u[1]))
    
    y = [u[1][0] - alpha]
    for j in range(1,N):
        y.append(y[j-1] + u[1][j] - (1+alpha)*u[0][j-1])

    # print(y)
    return y,ops

    
def durbins_split_linear_span_prog(N,p,k,u,U,U_dyn,ops,levinsons=False):
    # print("prog",p,k,u,U)
    # print(p[0],p[k],p[k-1])
    gamma = (p[0] + U_dyn[1][0] + p[k]) / (p[0] + U_dyn[0][0] + p[k-1])
    
    U_dyn,ops = durbins_split_linear_span_U_dyn_prog(N,p,k,gamma,U_dyn,ops,levinsons=levinsons)
        
    u_new_len = (k-1)//2
    # print(range(1,(k-1)//2+1))
    u_new = [u[1][0]-gamma+1]+[u[1][i] + u[1][i-1] - gamma*u[0][i-1] for i in range(1,(k-1)//2)]
    # for odd, looking back at last element misses end of half-list (handle as separate case)
    if k%2==1:
        u_new.append(u[1][-1] + u[1][-1] - gamma*u[0][-1])
    else:
        u_new.append(u[1][-1] + u[1][-2] - gamma*u[0][-1])
        
    # print("u_new",u_new)
    if k%2==0:
        # print(p[1:k//2+1],p[k//2+1:k+1])
        # print(u_new,rev(u_new))
        U_new = dot(p[1:k//2+1],u_new) + dot(p[k//2+1:k+1],rev(u_new))
    else:
        # print(p[1:k//2+2],p[k//2+2:k+1])
        # print(u_new,rev(u_new[:-1]))
        U_new = dot(p[1:k//2+2],u_new) + dot(p[k//2+2:k+1],rev(u_new[:-1]))

    # if k<N:
    #     print("U diff",U_dyn[1][0]-U_new)
    # print("U_new",U_new)
    # print("U_dyn new",U_dyn[1][0])
    ops=(ops[0]+k,ops[1]+k//2)
    return (u[1],u_new),(U[1],U_new),U_dyn,ops

def durbins_split_linear_span_U_dyn_prog(N,p,k,gamma,U_dyn,ops,levinsons=False):
    new_max=N-k if levinsons else N-k+1
    
    return (U_dyn[1],[U_dyn[1][q] + U_dyn[1][q+1] + p[1+q] + p[k+q] - (p[1+q]+p[k+q] + U_dyn[0][q+1])*gamma for q in range(0,new_max)]),(ops[0]+5*new_max,ops[1]+new_max)

                                                             
# u=[(y[0]+y[-1])/t[0]]
#     for k in range(2,N+1):
        
def durbins_split_linear_span(p,N):
    # print("p",p)
    # dubrin init
    y0=[-p[1]/p[0]]
    b0=p[0] + p[1]*y0[0]

    y1,b1,_ = durbins_prog(p,2,y0,b0,(0,0))
    y=(y0,y1)

    y_check=y1
    beta_check=b1
    
    # print(y)
    # we store all of u to make U initialization easier. Later, we only store first half of u (u is symmetric)
    u=tuple([y[i][j]+y[i][i-j] for j in range(i+1)] for i in range(2))
    # print("u",u)
    U=tuple(dot(u[i],p[1:i+2]) for i in range(2))
    # print("U",U)
 
    U_dyn = ([u[0][0]*p[q] for q in range(1,N+1)],[u[1][0]*p[q]+u[1][1]*p[q+1] for q in range(1,N)])
    # print("U_dyn",U_dyn)
    
    #debug
    

    # fix u so only half is stored
    u=(u[0],[u[1][0]])
    ops=(0,0)
    for k in range(3,N+1):
        # print("k",k)
        u,U,U_dyn,ops = durbins_split_linear_span_prog(N,p,k,u,U,U_dyn,ops)

        
    # recover y
    sum_from_half = lambda n,l : 2*sum(l) if n%2==0 else 2*sum(l)-l[-1]
    sum_u_n = sum_from_half(N,u[1])
    sum_u_nmo = sum_from_half(N-1,u[0])
    alpha = (sum_u_n - sum_u_nmo) / (sum_u_nmo + 2)

    u=(fill_from_half(N-1,u[0]),fill_from_half(N,u[1]))
    
    y = [u[1][0] - alpha]
    for j in range(1,N):
        y.append(y[j-1] + u[1][j] - (1+alpha)*u[0][j-1])

    # print(y)
    return y,ops
    

#

# solver :: n,row -> x,y
# durbins_generic :: t,y,N -> x

def durbins_generic_solver(durbins_generic):
    def solver(n,row):
        last = 100#2**(-n)
        y = -np.array(list(row[1:])+[last])
        x,ops = durbins_generic(list(row)+[last],n)
        # print(x)
        return x,y
    return solver


def durbins_generic_op_counter(durbins_generic):
    def counter(n,row):
        last = 100#2**(-n)
        y = -np.array(list(row[1:])+[last])
        x,ops = durbins_generic(list(row)+[last],n)
        # print(x)
        return ops
    return counter



    
if __name__=="__main__":
    import sys
    # dss = get_stability_sample(int(sys.argv[1]),10,kms,durbins_linear_span_solver)
    # print(dss)

    generics = [("durbins",durbins),
                ("durbins linear span",durbins_linear_span),
                ("durbins split",durbins_split),
                ("durbins split linear span",durbins_split_linear_span)]
                
    solvers = list(map(lambda ng : (ng[0],durbins_generic_solver(ng[1])),generics))
    counters = list(map(lambda ng : (ng[0],durbins_generic_op_counter(ng[1])),generics))
    

    # test correct
    
    n=int(sys.argv[1])
    sample = generate_sample(easy_gen,n,1)
    # determinants
    print(sample)
    print("determinants:\n   ",list(map(lambda m:np.linalg.det(m),sample)))
    
    for name,solver in [solvers[1]]:
        stab = get_stability(n,sample,solver)
        print(name,"\n   ",stab)


    
    # test numerical stability
    # print("### KMS ####")
    # stability_analysis(int(sys.argv[1]),100,kms,solvers)

    # print("### RND ####")
    # stability_analysis(int(sys.argv[1]),100,rnd,solvers)

    # print("### CVL ####")
    # stability_analysis(int(sys.argv[1]),100,cvl,solvers)

    
    # count ops

    n=int(sys.argv[1])
    sample = generate_sample(easy_gen,n,1)

    for s in sample:
        for name,counter in counters:
            print(name)
        
            count=counter(n,s[0])
            print("  add/sub:",count[0]/n**2)
            print("  mult:   ",count[1]/n**2)
