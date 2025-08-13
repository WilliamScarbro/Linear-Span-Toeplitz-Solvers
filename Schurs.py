import sys
from Util import *
from StabilityUtil import *

from Levinsons import levinsons_generic_solver,levinsons_generic_op_counter


def schur_prog(p,b,N,k,v,Vm,Vp,x,X):
    print("------------")
    phi = dot(p[1:k],v)
    print("phi ip",phi)

    Vp_ip = [dot(p[1+q:k+q],v) for q in range(N-k+1)]
    Vm_ip = [dot(rev(p[1+q:k+q]),v) for q in range(N-k+1)]

    # print("Vp",Vp)
    # print("Vp ip",Vp_ip)
    # print("Vm",Vm)
    # print("Vm ip",Vm_ip)
    # print("phi no ip",Vp[0])

    phi_ls = Vp[0]
    
    v_new_unscaled = [-phi_ls*v[k-2]] + [v[i-1]-phi_ls*v[k-2-i] for i in range(1,k-1)] + [v[k-2]]
    scale = 1/(1-phi**2)
    v_new = [scale*vnu for vnu in v_new_unscaled]

    gamma = Vp[0] # is this 1 or zero? math says 1, example says 0
    delta = 1/(1-gamma**2)
    # these appear to be flipped? 
    # Vp_new = [delta*(Vm[q+1]-gamma*Vp[q]) for q in range(N-k)]
    # Vm_new = [delta*(Vp[q]-gamma*Vm[q+1]) for q in range(N-k)]

    Vm_new = [delta*(Vp[q+1]-gamma*Vm[q]) for q in range(N-k)]
    Vp_new = [delta*(Vm[q]-gamma*Vp[q+1]) for q in range(N-k)]

    # assert T_k * v_new = e_0
    print("T_k v_k:",vec2toeplitz(p[0:k])@np.array(v_new))

    X_ip = dot(rev(p[1:k]),x)
    print("X_ip[0]",X_ip)
    eta = b[k-1] - X_ip
    # print("eta ip",eta)
    x_new = [x[i]+eta*v_new[i] for i in range(k-1)] + [eta*v_new[k-1]]

    X_ip = [dot(rev(p[1+q:k+q]),x) for q in range(N-k+1)]
    print("X ip",X_ip)
    print("X",X)

    eta_no_ip = b[k-1] - X[0]
    # print("eta no ip",eta_no_ip)
    X_new = [X[q] + eta_no_ip*Vm[q] for q in range(N-k)]

    # assert T_k * x_new = b_k
    # print("T_k x_k:",vec2toeplitz(p[0:k])@np.array(x_new))
    # print("b_k",b[:k])

    
    return v_new,Vp_new,Vm_new,x_new,X_new


def schur_levinsons(p,b,N):
    print("p",p)
    print("b",b)
    
    # v init (T_1 v_1 = 1)
    v=[1/p[0]]
    
    # x init (T_1 x_1 = b_0)
    x=[b[0]/p[0]]

    # V init
    Vp = [v[0]*p[i] for i in range(1,N)]
    Vm = Vp.copy()
    
    # X init
    X = [x[0]*p[i] for i in range(1,N)]
    
    for k in range(2,N+1):
        v,Vp,Vm,x,X = schur_prog(p,b,N,k,v,Vp,Vm,x,X)

    # probably requires last step on only x
    
    return x,0

if __name__=="__main__":
    sl_solver = levinsons_generic_solver(schur_levinsons)
    sl_counter = levinsons_generic_op_counter(schur_levinsons)

    n = int(sys.argv[1])
    sample = generate_sample(easy_gen,n,1)
    print(sample)
    stab = get_stability(n,sample,sl_solver)
    print("error:",stab)

    
