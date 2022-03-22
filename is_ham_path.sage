def is_ham_path(bcount, m, n, e):
    r"""In the abelian Cayley digraph `Cay(G; a, b)`, such that `o(a) = m`, 
    `|G : <a>| = n`, and `nb = ea`, construct a spanning quasipath `H` in which
    the number of vertices that travel by `b` is `bcount`.  Return `True` if `H`
    is a hamiltonian path."""

    def E(v):
        r"""Return the vertex `v + (1,0)`."""
        return ( (v[0] + 1) % m, v[1])
    def W(v):
        r"""Return the vertex `v - (1,0)`."""
        return ( (v[0] - 1) % m, v[1])
    def N(v):
        r"""Return the vertex `v + (0,1)`."""
        if v[1] != n - 1:
            return (v[0], v[1] + 1)
        else:
            return ((v[0] + e) % m, 0)
    def S(v):
        r"""Return the vertex `v - (0,1)`."""
        if v[0] != 0:
            return (v[0], v[1] - 1)
        else:
            return ((v[0] - e) % m, n - 1)
    def NW(v):
        r"""Return `N(W(v))`.  We use this function to iterate through the
        arc-forcing coset that contains `v`."""
        return N(W(v))
    
    number_of_cosets = gcd(m, n - e)  # number of cosets of the arc-forcing subgroup
    arc_forcing_order = m * n // number_of_cosets  # order of the arc-forcing subgroup
    t, d = bcount.quo_rem(arc_forcing_order)
    init = (0,0)  # the initial vertex of `H`

    # We represent the spanning quasi-path `H_t(d)` as a dictionary: `H[v]` is
    # the unique out-neighbour of the vertex `v` (except that the terminal vertex
    # has no out-neighbour). We start with an empty dictionary, and add the
    # dictionary entries one-by-one.

    H = {}

    # the first `t` non-terminal cosets travel by `(0,1)`
    # and the rest travel by `(1,0)`
    coset_rep = init
    for i in range(number_of_cosets - 1):
        v = coset_rep
        for _ in range(arc_forcing_order):
            H[v] = N(v) if i < t else E(v)
            v = NW(v)
        coset_rep = E(coset_rep)
    
    # the first `d` vertices in the terminal coset travel by (0,1)
    v = W(init)
    for _ in range(d):
        H[v] = N(v)
        v = NW(v)
    # the next vertex is the terminal vertex (which has no out-neighbour)
    term = v
    # the rest of the vertices in the terminal coset travel by (1,0)
    for _ in range(arc_forcing_order - d - 1):
        v = NW(v)
        H[v] = E(v)    

    # To determine whether `H` is a hamiltonian path, we check whether the path
    # component of `H` has the correct length    
    path_length = 0
    v = init
    while v != term:
        path_length += 1
        v = H[v]
    return path_length == m * n - 1
