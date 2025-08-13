def dot(a,b):
    assert len(a)==len(b), f"{a},{b}"
    return sum(ai*bi for ai,bi in zip(a,b))


def rev(a):
    b=list(a.copy())
    b.reverse()
    return b

def fill_from_half(k,half):
    if k%2==0:
        return half+rev(half)
    return half+rev(half[:-1])

