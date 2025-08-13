
import math
import numpy as np
import random

def vec2toeplitz(x,symmetric=True):
    return np.array([[x[abs(i-j)] for i in range(len(x))] for j in range(len(x))])

def cvl(theta,n):
    t = lambda theta : [math.cos(2*math.pi*theta*i) for i in range(n)]
    thetas = [random.random() for i in range(n)]
    ws = [random.random() for i in range(n)]

    T = [sum(ws[k]*t(thetas[k])[i] for k in range(n)) for i in range(n)]
    m = 1/T[0]
    mT = [m*T[i] for i in range(n)]
    return vec2toeplitz(mT)


def kms(nu,n):
    x=[nu**abs(i) for i in range(n)]
    return vec2toeplitz(x)

def rnd(_,n):
    x=[random.random() for i in range(n)]
    return vec2toeplitz(x)

def easy_gen(_,n):
    #x=[1.0*2**(-i) for i in range(n)]
    x=[i+1 for i in range(n)]
    return vec2toeplitz(x)

def generate_sample(gen_func,mat_size,sample_size):
    return [gen_func(random.random(),mat_size) for i in range(sample_size)]

def l2(arr):
    return np.power(np.sum(np.power(arr,2)),1/2)


def get_stability(n,sample,solver):
    def solve(s):
        return s,*solver(n,s[0])
    def get_error(s,x,y):
        # print("expected",y)
        # print("actual",s@x)
        # print("diff",s@x-y)
        true=s@x
        return l2(true - y)/l2(true)
    
    return list(map(lambda s : get_error(*solve(s)),sample))

def get_stability_sample(n,sample_size,gen_func,solver):
    sample = generate_sample(gen_func,n,sample_size)
    return get_stability(n,sample,solver)


def stability_analysis(N,size,generator,solvers):
    sample=generate_sample(generator,N,size)


    def average(l):
        return sum(l)/len(l)

    def log(l):
        return np.array(list(map(lambda x:math.log(x+1e-24)/math.log(10),l)))

    means = []
    for name,solver in solvers:
        errors = get_stability(N,sample,solver)
        # print(name)
        # print("   mean: ",average(errors))
        # print("   geo mean: ",10**np.mean(log(errors)))
        # print("   geo std:  ",10**np.std(log(errors)))
        print(np.mean(log(errors)),",",np.std(log(errors)))
        # means.append(np.mean(log(errors)))
        # print("   diff:     ",means[0]-means[-1])

    # mean_diffs = [mean[0]-mean[i] for i in range(1,len(solvers))]
    

        
def test_solver(n,row):
    return np.zeros(n),row


if __name__=="__main__":
    get_stability_sample(4,10,kms,test_solver)

    

