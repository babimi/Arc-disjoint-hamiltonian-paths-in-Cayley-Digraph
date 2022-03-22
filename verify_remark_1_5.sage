# Running the sagemath code in this file will verify Remark 1.5 of the paper
# "Arc-disjoint hamiltonian paths in Cartesian products of directed cycles"
# by Iren Darijani, Babak Miraftab, and Dave Witte Morris.  It establishes
# that if G is an abelian group of order at most 10^4, then every 2-generated,
# connected, loopless Cayley digraph on G has two arc-disjoint hamiltonian paths.
# By Proposition 6.2 of the paper, only two specific families of Cayley digraphs
# need to be considered.

# The following file does not need to be loaded if the `confirm_paths` argument
# is set to `False` (which is the default) in the functions `has_arc_disjoint`
# and `remark_1_5`.
load("is_ham_path.sage")

def prim_latt_pts(m, n, e):
    r"""Return a list of the primitive lattice points in the closed triangle
    T(m, n, e), except the first one (on the x-axis) and the last one (on the
    side from (0, 0) to (e, m). They are sorted in order of increasing slope."""
    result = []
    for y in range(1, m):
        for x in range(((y*e) // m) + 1, ((m*n - (n - e)*y) // m) + 1):
            if (Integer(x).gcd(y) == 1):
                result.append((x,y))
    result.sort(key=lambda p: p[1]/p[0])
    return result

def bcount_list(m, n, e):
    r"""For each hamiltonian path in the abelian Cayley digraph Cay(G; a, b)
    with parameters m, n, e, record the number of vertices that travel by b in
    the hamiltonian path. Return these numbers in a sorted list (from smallest
    to largest)."""
    result = [n-1]
    for x,y in prim_latt_pts(m, n, e):
        delta = 2 * ((m*n - 1)//(m*x + (n - e)*y))
        if (m*n) % (m*x + (n - e)*y) == 0:
            delta += 1
        result.append(result[-1] + delta)
    return result

def has_arc_disjoint(m, n, e, confirm_paths=False):
    r"""Return `True` if the abelian Cayley digraph Cay(G; a, b) with
    parameters m, n, e has two arc-disjoint hamiltonian paths.

    Setting `confirm_paths=True` will confirm that the hamiltonian paths actually
    exist, instead of relying on Theorem 3.21 (the Curran-Witte Theorem), but
    takes more than twice as long."""
    numb_list = bcount_list(m, n, e)
    for numb in numb_list:
        opposite = m*n - 1 - numb
        if ((opposite in numb_list) or (opposite - 1 in numb_list)
            or (opposite + 1 in numb_list)):
                if confirm_paths:
                    numb_opp = (opposite if opposite in numb_list
                                else opposite - 1 if opposite - 1 in numb_list
                                else opposite + 1)
                    assert (is_ham_path(numb, m, n, e) 
                            and is_ham_path(numb_opp, m, n, e))
                return True
    return False

def remark_1_5(min_order, max_order, show_progress=0, confirm_paths=False):
    r"""Consider all Cayley digraphs listed in Proposition 6.2 whose order is
    between `min_order` and `max_order` (inclusive). Print a list of the values
    of `m`, `n`, `e` for which there does not exist a pair of arc-disjoint
    hamiltonian paths.

    Setting the keyword `show_progress` to a nonzero integer value `x` will
    cause a message to be printed whenever the order is a multiple of `x`."""
    bad_ones = []
    if show_progress:
        from time import time
        starttime = time()
    for k in srange(min_order, max_order + 1):
        for a in srange(1, k // 2 + 1):

            # part (1) of Proposition 6.2
            b = a + 1
            n = a.gcd(k)
            m = k // n
            e = (b * (a // n).inverse_mod(m)) % k
            assert (n*b - e*a) % k == 0
            if is_odd(gcd(m, n - e)) and not(has_arc_disjoint(m, n, e, confirm_paths)):
                bad_ones.append([m, n, e])
                print(f"    no arc-disjoint hamiltonian paths for ({m}, {n}, {e})")

            # part (2) of Proposition 6.2
            if k % (2*a + 1) == 0:
                b = a + 1
                n = a.gcd(k)
                m = k // n
                e = (b * (-a//n).inverse_mod(m)) % k
                assert (n*b - e*(-a)) % k == 0
                if is_odd(gcd(m, n - e)) and not(has_arc_disjoint(m, n, e, confirm_paths)):
                    bad_ones.append((m, n, e))
                    print(f"    no arc-disjoint hamiltonian paths for ({m}, {n}, {e})")
                
        if show_progress and (k % show_progress == 0): 
            print(f"finished order {k}: elapsed time = {round((time() - starttime)/60, 2)} minutes")

    if bad_ones:
        print(f"\nSome (m, n, e) do not have arc-disjoint hamiltonian paths:\n{badones}")
    else:
        print("\nDone. All have arc-disjoint hamiltonian paths.")

# Verify Remark 1.5:
remark_1_5(3, 10^4, show_progress=10, confirm_paths=False)
