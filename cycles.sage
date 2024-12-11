
# Functions for navigating isogeny graphs and finding cycles.
from sage.schemes.elliptic_curves.ell_curve_isogeny import compute_isogeny_bmss

def get_neighbors(j, ell):
    Phi_inst = classical_modular_polynomial(ell, j)
    neighbors = Phi_inst.roots()
    return neighbors

def adjacent_vertices(j, ell):
    """
    Returns the j-invariants adjacent to j in the ell-isogeny graph.
    """
    return [root for root,mult in get_neighbors(j, ell) for _ in range(mult)]

def random_walk(j0, ell):
    """
    Returns a random walk in the ell-isogeny graph rooted at j0.
    
    j0 is assumed to be then j-invariant of an elliptic curve over a finite field of characteristic p containing Fp2. In this case, the walk of length log(p) will terminate in an approximatley uniformly distributed supersingular j-invariant."""
    p = j0.parent().characteristic()
    assert p > 0
    length = floor(log(p, ell))
    path = [j0]
    for k in range(length):
        jk = path[k]
        neighbors = adjacent_vertices(jk, ell)
        path.append(choice(neighbors))
    return path

def random_walks_to_Fp(j0, ell):
    """Returns a non-backtracking path in the ell-isogeny graph from j0 to a random supersingular j-invariant in Fp."""
    p = j0.parent().characteristic()
    length = floor(log(p, ell))
    while True:
        walk = [randint(0, ell)] + [randint(0, ell-1) for _ in range(length-1)]
        j1 = adjacent_vertices(j0, ell)[walk[0]]
        path = [j0, j1]
        for k in range(1,length):
            jk = path[k]
            neighbors = adjacent_vertices(jk, ell)
            # remove jk from neighbors so we don't backtrack. jk may still be in neighbors if there are multiple edges at jk.
            neighbors.remove(path[k-1])
            j = neighbors[walk[k]]
            path.append(j)
        if path[-1]**p == path[-1]:
            return path


def reflect_path(path):
    """Returns the image of the edges in path under the Frobenius isomorphism on the isogeny graph G(p,ell)"""
    kernel_polys = [phi.kernel_polynomial() for phi in path]
    R = kernel_polys[0].parent()
    p = path[0].domain().base_field().characteristic()
    galois_conjugate_polys = [R([a**p for a in h.coefficients()]) for h in kernel_polys]
    galois_conjugate_curves = [frob(phi.domain()) for phi in path]
    return [EllipticCurveIsogeny(galois_conjugate_curves[k], galois_conjugate_polys[k]) for k in range(len(path))]

def CGL_cycle(j, ell):
    """ Compute a cycle in G(p,ell) at j.

    When j is in Fp, this cycle may be trivial.
    """
    p = j.parent().characteristic()
    P1 = random_walks_to_Fp(j, ell)
    P1_conjugate = [j**p for j in P1]
    P1_conjugate.reverse()
    half_cycle_1 = P1 + P1_conjugate[1:]
    #if j in Fp:
    #   return half_cycle_1
    P2 = random_walks_to_Fp(j, ell)
    while P2[-1] == P1[-1]:
        P2 = random_walks_to_Fp(j, ell)
    P2_conjugate = [j**p for j in P2]
    P2.reverse()
    half_cycle_2 = P2_conjugate[1:] + P2[1:]
    cycle = half_cycle_1 + half_cycle_2
    return cycle


def path_to_isogeny(path, ell):
    """
    Given path a list of adjacent j-invariants in G(p,ell),
    return a corresponding (factored) isogeny.

    If the codomain of the chain is isomorphic to the domain,
    an endomorphism is returned.
    """
    if len(path) < 2:
        raise NotImplementedError
    E = EllipticCurve(j=path[0])
    chain = []
    for j in path[1:]:
        Anext, Bnext = normalized_model(E,j,ell)
        Enext = EllipticCurve([Anext, Bnext])
        phi = E.isogeny(None, Enext, ell) if ell == 2 else E.isogeny(compute_isogeny_bmss(E, Enext, ell))
        chain.append(phi)
        E = chain[-1].codomain()
    try:
        iso = E.isomorphism_to(chain[0].domain())
    except ValueError:
        pass
    else:
        chain[-1] = iso * chain[-1]
    from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite   # great opportunity to update Sage!
    return EllipticCurveHom_composite.from_factors(chain)

def find_supersingular_j(p):
    """ Returns a supersingular j-invariant in Fp2 using Broker's algorithm."""
    if not p.is_prime():
        raise ValueError('input must be prime')
    F = GF(p^2)
    if p%12 == 7:
        # E = EllipticCurve([F(1),F(0)])
        j = F(1728)
    elif p%12 == 5:
        # E = EllipticCurve([F(0),F(1)])
        j = F(0)
    elif p%12 == 11:
        # E = EllipticCurve([F(1),F(0)])
        j = F(1728)
    else:
        q = 3
        while kronecker(-q,p) == 1 or q%4 == 1:
            q = next_prime(q)
        PK = hilbert_class_polynomial(-q)
        Fx = F['x']
        PK = Fx(PK)
        j = PK.roots()[0][0]
    return j
