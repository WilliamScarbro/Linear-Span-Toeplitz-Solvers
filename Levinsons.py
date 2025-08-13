import Durbins
from StabilityUtil import *
from Util import *

# assumes symmetric T (t in R^{n})

def levinsons_prog(p,b,k,y,x,ops):
    # print('lev prog')
    # print("p",p)
    # print("b",b)
    # print("y",y)
    # print("x",x)
    mu = (b[k-1] - dot(p[1:k],rev(x))) / (p[0] + dot(p[1:k],y))
    # print("mu",mu)
    x_new = [x[i] + mu*y[k-2-i] for i in range(k-1)]+[mu]
    # print("x_new&mu check",dot(x_new[:-1],rev(p[1:k]))+mu*p[0],b[k-1])
    ops = (ops[0]+k+k,ops[1]+k+k)
    return x_new,ops

def levinsons(p,b,N):
    x=[b[0]/p[0]]
    y=[-p[1]/p[0]]
    d=p[0] + p[1]*y[0] # beta in durbin's algorithm
    # p=list(p)+[0]

    # print(vec2toeplitz(p[0:1])@np.array(y)+p[1])
    ops=(0,0)
    for k in range(2,N):
        # print(k,len(x),len(y))
        x,ops   = levinsons_prog(p,b,k,y,x,ops)
        y,d,ops = Durbins.durbins_prog(p,k,y,d,ops)
        # print("y check",vec2toeplitz(p[0:k])@np.array(y)+np.array(p[1:k]))

        # print(k,len(x),len(y))

    x,ops = levinsons_prog(p,b,N,y,x,ops)
    return x,ops



def levinsons_linear_span_X_prog(N,p,b,k,mu,X_hat,Y,ops):
    return [X_hat[q+1]+mu*(Y[q+1]+p[1+q]) for q in range(0,N-k-1)],(ops[0]+2*(N-k),ops[1]+N-k)

def levinsons_linear_span_prog(N,p,b,k,y,x,X_hat,Y,ops):
    # print("X_hat",X_hat)
    # print("X_hat true",[dot(p[1+q:k+q],rev(x)) for q in range(N-k+1)])
    # print("Y",Y[0])
    # print("Y true",dot(p[1:k],y))
    # mu_true = (b[k-1] - dot(p[1:k],rev(x))) / (p[0] + dot(p[1:k],y))
    mu = (b[k-1] - X_hat[0]) / (p[0] + Y[0])
    # print("mu",mu)
    # print("mu_true",mu_true)
    
    X_hat_new,ops = levinsons_linear_span_X_prog(N,p,b,k,mu,X_hat,Y,ops)
    
    # print("mu",mu)
    x_new = [x[i] + mu*y[k-2-i] for i in range(k-1)]+[mu]
    # print("x_new&mu check",dot(x_new[:-1],rev(p[1:k]))+mu*p[0],b[k-1])
    ops=(ops[0]+k,ops[1]+k)
    return x_new,X_hat_new,ops
    
def levinsons_linear_span(p,b,N):
    #durbin init
    y=[-p[1]/p[0]]
    Y_hat=[p[1+q]*y[0] for q in range(N-1)]
    Y=Y_hat.copy()
   
    # levinson init
    x=[b[0]/p[0]]
    X_hat=[x[0]*p[1+q] for q in range(N-1)]
    # print("X_hat",X_hat)
    

    ops=(0,0)
    for k in range(2,N):
        x,X_hat,ops = levinsons_linear_span_prog(N,p,b,k,y,x,X_hat,Y,ops)
        Y,Y_hat,y,ops = Durbins.durbins_linear_span_prog(N,p,k,Y,Y_hat,y,ops,for_levinsons=True)

    x,ops = levinsons_prog(p,b,N,y,x,ops)
    return x,ops
        

def levinsons_split_prog(N,p,b,k,w,u,W,U,ops):
    gamma = (W[1] - W[0] - b[k-1] + b[k-2]) / (p[0] + U[0] + p[k-1])

    w_new = [w[1][0] - gamma] + [w[1][i] + w[1][i-1] - w[0][i-1] - gamma*u[0][i-1] for i in range(1,(k-1)//2)]
    W_new=0
    if k%2==1:
        w_new.append(w[1][-1] + w[1][-1] - w[0][-1] - gamma*u[0][-1])
        # print("w_new",w_new)
        if k < N:
            # print(p[1:k//2+2],p[k//2+2:k+1])
            # print(w_new,rev(w_new[:-1]))
            W_new = dot(p[1:k//2+2],w_new) + dot(p[k//2+2:k+1],rev(w_new[:-1]))
        
    else:
        w_new.append(w[1][-1] + w[1][-2] - w[0][-1] - gamma*u[0][-1])
        # print("w_new",w_new)
        if k<N:
            # print(p[1:k//2+1],p[k//2+1:k+1])
            # print(w_new,rev(w_new))
 
            W_new = dot(p[1:k//2+1],w_new) + dot(p[k//2+1:k+1],rev(w_new))

    ops=(ops[0]+3*(k//2)+k,ops[1]+k//2+k//2)
    return (w[1],w_new),(W[1],W_new),ops

def levinsons_split(p,b,N):
    ops=(0,0)
    # print("p",p)
    # print("b",b)
    
    # durbins init
    y0=[-p[1]/p[0]]
    b0=p[0] + p[1]*y0[0]

    y1,b1,ops = Durbins.durbins_prog(p,2,y0,b0,ops)
    y=(y0,y1)

    y_check=y1
    beta_check=b1
    
    # print(y)
    # we store all of u to make U initialization easier. Later, we only store first half of u (u is symmetric)
    u=tuple([y[i][j]+y[i][i-j] for j in range(i+1)] for i in range(2))
    # print("u",u)
    U=tuple(dot(u[i],p[1:i+2]) for i in range(2))
    # print("U",U)


    # levinsons init
    x0 = [b[0]/p[0]]
    x1,ops = levinsons_prog(p,b,2,y0,x0,ops)
    # print("x",(x0,x1))

    w = ([2*x0[0]],[x1[0]+x1[1]]*2)
    # print("w",w)
    W = (dot(w[0],p[1:2]), dot(w[1],p[1:3]))
    # print("W",W)

    u=(u[0],[u[1][0]])
    w=(w[0],[w[1][0]])

    for k in range(3,N):
        # print("-------\nk",k)
        w,W,ops = levinsons_split_prog(N,p,b,k,w,u,W,U,ops)
        u,U,ops2 = Durbins.durbins_split_prog(p,k,u,U,ops)

    w,W,ops = levinsons_split_prog(N,p,b,N,w,u,W,U,ops)
    # last u cannot be computed (runs out of p)
    u=(u[1],None)
    
    # recover x
    sum_from_half = lambda n,l : 2*sum(l) if n%2==0 else 2*sum(l)-l[-1]
    sum_w_n = sum_from_half(N,w[1])
    sum_w_nmo = sum_from_half(N-1,w[0])
    sum_u_nmo = sum_from_half(N-1,u[0])
    mu = (sum_w_n - sum_w_nmo) / (sum_u_nmo + 2)

    w=(fill_from_half(N-1,w[0]),fill_from_half(N,w[1]))
    u=(fill_from_half(N-1,u[0]),None)

    # print("w full",w)
    # print("u full",u)
    
    x = [w[1][0] - mu]
    for j in range(1,N):
        x.append(x[j-1] + w[1][j] - w[0][j-1] - mu*u[0][j-1])

    return x,ops




######  split linear span

def levinsons_split_linear_span_W_dyn_prog(N,p,k,gamma,U,W,ops):
    return (W[1],[W[1][q] + W[1][q+1] - W[0][q+1] - (U[0][q+1]+p[1+q]+p[k+q])*gamma for q in range(0,N-k)]),(ops[0]+5*(N-k),ops[1]+N-k)

def levinsons_split_linear_span_prog(N,p,b,k,w,u,W,U,U_dyn,W_dyn,ops):
    gamma = (W_dyn[1][0] - W_dyn[0][0] - b[k-1] + b[k-2]) / (p[0] + U_dyn[0][0] + p[k-1])

    w_new = [w[1][0] - gamma] + [w[1][i] + w[1][i-1] - w[0][i-1] - gamma*u[0][i-1] for i in range(1,(k-1)//2)]
    W_dyn_new,ops = levinsons_split_linear_span_W_dyn_prog(N,p,k,gamma,U_dyn,W_dyn,ops)
    
    W_new=0
    if k%2==1:
        w_new.append(w[1][-1] + w[1][-1] - w[0][-1] - gamma*u[0][-1])
        # print("w_new",w_new)
        if k < N:
            # print(p[1:k//2+2],p[k//2+2:k+1])
            # print(w_new,rev(w_new[:-1]))
            W_new = dot(p[1:k//2+2],w_new) + dot(p[k//2+2:k+1],rev(w_new[:-1]))
        
    else:
        w_new.append(w[1][-1] + w[1][-2] - w[0][-1] - gamma*u[0][-1])
        # print("w_new",w_new)
        if k<N:
            # print(p[1:k//2+1],p[k//2+1:k+1])
            # print(w_new,rev(w_new))
 
            W_new = dot(p[1:k//2+1],w_new) + dot(p[k//2+1:k+1],rev(w_new))

    # if k<N:
    #     # print("W",W_new)
    #     print("W-W_dyn",W_dyn[1][0]-W_new)
    ops=(ops[0]+3*(k//2),ops[1]+k//2)
    return (w[1],w_new),(W[1],W_new),W_dyn_new,ops
    
def levinsons_split_linear_span(p,b,N):
    ops=(0,0)
    # print("p",p)
    # print("b",b)
    
    # durbins init
    y0=[-p[1]/p[0]]
    b0=p[0] + p[1]*y0[0]

    y1,b1,ops = Durbins.durbins_prog(p,2,y0,b0,ops)
    y=(y0,y1)

    y_check=y1
    beta_check=b1
    
    # print(y)
    # we store all of u to make U initialization easier. Later, we only store first half of u (u is symmetric)
    u=tuple([y[i][j]+y[i][i-j] for j in range(i+1)] for i in range(2))
    # print("u",u)
    U=tuple(dot(u[i],p[1:i+2]) for i in range(2))
    # print("U",U)
    U_dyn = ([u[0][0]*p[q] for q in range(1,N)],[u[1][0]*p[q]+u[1][1]*p[q+1] for q in range(1,N-1)])

    
    # levinsons init
    x0 = [b[0]/p[0]]
    x1,ops = levinsons_prog(p,b,2,y0,x0,ops)
    # print("x",(x0,x1))

    w = ([2*x0[0]],[x1[0]+x1[1]]*2)
    # print("w",w)
    W = (dot(w[0],p[1:2]), dot(w[1],p[1:3]))
    # print("W",W)

    W_dyn = ([w[0][0]*p[q] for q in range(1,N)],[w[1][0]*p[q]+w[1][1]*p[q+1] for q in range(1,N-1)])
    # print("W_dyn",W_dyn)
    
    u=(u[0],[u[1][0]])
    w=(w[0],[w[1][0]])

           
    for k in range(3,N):
        # print("-------\nk",k)
        w,W,W_dyn,ops = levinsons_split_linear_span_prog(N,p,b,k,w,u,W,U,U_dyn,W_dyn,ops)
        u,U,U_dyn,ops = Durbins.durbins_split_linear_span_prog(N,p,k,u,U,U_dyn,ops,levinsons=True)

    w,W,W_dyn,ops = levinsons_split_linear_span_prog(N,p,b,N,w,u,U,W,U_dyn,W_dyn,ops)
    
    # last u cannot be computed (runs out of p)
    u=(u[1],None)
    
    # recover x
    sum_from_half = lambda n,l : 2*sum(l) if n%2==0 else 2*sum(l)-l[-1]
    sum_w_n = sum_from_half(N,w[1])
    sum_w_nmo = sum_from_half(N-1,w[0])
    sum_u_nmo = sum_from_half(N-1,u[0])
    mu = (sum_w_n - sum_w_nmo) / (sum_u_nmo + 2)

    w=(fill_from_half(N-1,w[0]),fill_from_half(N,w[1]))
    u=(fill_from_half(N-1,u[0]),None)

    # print("w full",w)
    # print("u full",u)
    
    x = [w[1][0] - mu]
    for j in range(1,N):
        x.append(x[j-1] + w[1][j] - w[0][j-1] - mu*u[0][j-1])

    return x,ops




#

# solver :: n,row -> x,y
# levinsons_generic :: t,y,N -> x

def levinsons_generic_solver(levinsons_generic):
    def solver(n,row):
        b = np.array([-2**(-i-1) for i in range(n)])
        x,ops = levinsons_generic(row,b,n)
        return x,b
    return solver

def levinsons_generic_op_counter(levinsons_generic):
    def counter(n,row):
        b = np.array([-2**(-i-1) for i in range(n)])
        x,ops = levinsons_generic(row,b,n)
        return ops
    return counter





if __name__=="__main__":
    
    import sys
    
    generics = [("levinsons",levinsons),
               ("levinsons linear span",levinsons_linear_span),
               ("levinsons split",levinsons_split),
               ("levinsons split linear span",levinsons_split_linear_span)]

    solvers = list(map(lambda ng : (ng[0],levinsons_generic_solver(ng[1])),generics))
    counters = list(map(lambda ng : (ng[0],levinsons_generic_op_counter(ng[1])),generics))

    # n=int(sys.argv[1])
    # sample = generate_sample(easy_gen,n,1)
    # # determinants
    # print(sample)
    # print("determinants:\n   ",list(map(lambda m:np.linalg.det(m),sample)))
    
    # for name,solver in solvers[0:4]:
    #     stab = get_stability(n,sample,solver)
    #     print(name,"\n   ",stab)
        
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
        for name,counter in counters[0:4]:
            print(name)
        
            count=counter(n,s[0])
            print("  add/sub:",count[0]/n**2)
            print("  mult:   ",count[1]/n**2)
