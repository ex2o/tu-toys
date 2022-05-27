N = 10

b_shapley = [0 if _n<(N-1) else 1 for _n in range(N)]
b_banzhaf = [1/(2**(N-1)) for _ in range(N)]

b = b_shapley

def b_to_weights(b):
    n = len(b)
    #w = [None for _ in range(n)]
    # Not python: w[n] = b[n]/n
    w = [_b/(_n+1) for _n, _b in enumerate(b)]
    for t in range(n-2, -1, -1):
        w[t] = (b[t] + (n-t-1)*w[t+1])/(t+1)

    return w

print(b_to_weights(b))
