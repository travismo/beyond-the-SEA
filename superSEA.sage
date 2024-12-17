def proj_double(P, A, f, h):
    """Doubles the projective point P with coordinates in F[x]/h(x).
     
      P: the restriction of an endomorphism alpha to the kernel of an isogeny with kernel polynomial h. P is representated by polynomials u(x), v(x), w(x) so that if 
      alpha(X,Y,Z) = [U(X,Z) : V(X,Z)Y: W(X,Z)Z] in projective coordinates, then 
      u(x) = U(x,1) mod h, v(x) = V(x,1) mod h, w(x) = W(x,1) mod h.
      A: a4 coefficient of domain of P
      f: domain of P is y^2=f(x)
      h: a kernel polynomial for an ell isogeny for some odd prime ell 
     """
    (x1, y1, z1) = P
    # x3 = 2 * y1 * z1 * ((A*z1^2 + 3*x1^2)^2 - 8 * x1 * y1^2 * z1)
    # y3 = (3 * x1^2  + A * z1^2) * (12 * x1 * y1^2 * z1 - (3 * x1^2 + A * z1^2)^2) - 8 * y1^4 * z1^2
    # z3 = (2 * y1 * z1)^3
    # multiply above through by y, use y^2z = f(x,z) so output has form (U(X,Z):V(X,Z)Y:W(X,Z)Z) with U,V,W not divisible by z
    fprime = 3 * x1^2 + A * z1^2 % h 
    fprimesq = fprime^2 % h 
    y1z1 = y1 * z1 % h 
    y1z1f = y1z1 * f % h 
    y1z1fsq = y1z1f^2 % h 
    x1y1sqz1f = x1 * y1  % h * y1z1f % h 
    x3 = 2 * y1z1f  % h * (fprimesq - 8 * x1y1sqz1f) % h 
    y3 = fprime * (12 * x1y1sqz1f - fprimesq) % h  - 8 * y1^2  % h * y1z1fsq % h 
    z3 = 8 * y1z1 * y1z1fsq % h 
    #### TODO reuse intermediate multiplications/squarings
    return (x3, y3, z3)
def proj_add(P, Q, A, f, h):
    """Adds the projective points P and Q with coordinates in F[x]/h(x).
     
      P: the restriction of an endomorphism alpha to the kernel of an isogeny with kernel polynomial h. P is representated by polynomials u(x), v(x), w(x) so that if 
      alpha(X,Y,Z) = [U(X,Z) : V(X,Z)Y: W(X,Z)Z] in projective coordinates, then 
      u(x) = U(x,1) mod h, v(x) = V(x,1) mod h, w(x) = W(x,1) mod h.
      Q: the restriction of an endomorphism to the kernel of an isogeny with kernel polynomial h.
      A: a4 coefficient of domain of P
      f: domain of P is y^2=f(x)
      h: a kernel polynomial for an ell isogeny for some odd prime ell 
     """
    if isequal(P, Q, h): return proj_double(P, A, f, h)
    (x1, y1, z1) = P
    (x2, y2, z2) = Q
    x2z1 = x2 * z1 % h 
    x1z2 = x1 * z2 % h 
    x2z1mx1z2 = x2z1 - x1z2  % h 
    x2z1mx1z2sq = x2z1mx1z2^2 % h 
    x2z1mx1z2cu = x2z1mx1z2sq * x2z1mx1z2 % h 
    x2z1px1z2 = x2z1 + x1z2 % h 
    z1z2 = z1 * z2 % h 
    z1z2f = z1z2 * f % h 
    y1z2 = y1 * z2 % h 
    y2z1 = y2 * z1 % h 
    y2z1my1z2 = y2z1 - y1z2 % h 
    y2z1my1z2sq = y2z1my1z2^2 % h 
    x3 = (x2z1mx1z2) * (y2z1my1z2sq * z1z2f  % h - (x2z1mx1z2sq) * (x2z1px1z2) % h ) % h 
    y3 = y2z1my1z2 * (x2z1mx1z2sq * (x2z1px1z2 + x1z2) % h - y2z1my1z2sq * z1z2f % h )  % h - (x2z1mx1z2cu) * y1z2 % h 
    z3 = (x2z1mx1z2cu) * z1z2 % h 
    ### TODO reuse intermediate multiplications/squarings
    return (x3, y3, z3)

def negate(P, h):
    (x, y, z) = P
    return (x, -y % h, z)

def isequal(P, Q, h):
    """Equality testing for restrictions of endomorphisms in projective coordinates."""
    (x1, y1, z1) = P
    (x2, y2, z2) = Q
    if x1 * y2 % h != x2 * y1 % h: return False
    if y1 * z2 % h != y2 * z1 % h: return False
    return True

def iszero(P, h):
    return isequal(P, [0,1,0], h)

def scale(n, P, A, f, h):
    """
    Compute the scalar multiple n*P in Hom(ker phi, E[ell]) using double and add.
    """
    if not n: return ()
    nbits = n.digits(2)
    i = len(nbits) - 2
    Q = P
    while i >= 0:
        Q = proj_double(Q, A, f, h)
        if nbits[i]: Q = proj_add(P, Q, A, f, h)
        i -= 1
    return Q

def projective_chain(endo):
    """
    Compute functions u(x), v(x), w(x) for each factor phi of endo so that if 
      phi(X,Y,Z) = [U(X,Z) : V(X,Z)Y: W(X,Z)Z] in projective coordinates, then 
      u(x) = U(x,1), s(x) = V(x,1), d(x) = W(x,1)
    """
    chain = []
    for phi in endo.factors():
        phix = phi.x_rational_map()
        u = phix.numerator()
        v = phix.denominator()
        s = u.derivative()*v - u*v.derivative()
        t = phi.scaling_factor() * v^2
        d = v.lcm(t)
        u = u * (d // v)
        s = s * (d // t)
        # phi(X:Y:Z) = [U(X,Z) : S(X,Z)Y : D(X,Z)Z], where U,S,D are the homogenizations of u(x),s(x),d(x). 
        chain.append((u,s,d))
    return chain
    
def proj_isogeny_compose_mod(phi, alpha, h):
    """Returns the coordinate functions of phi*alpha modulo h. 
    
    phi: a triple (u, s, d) of polynomials in x representing an isogeny in projective form given by [U(X,Z) : S(X,Z)Y : D(X,Z)Z] where U,S,D are the homogenizations of u,s,d
    alpha: a triple (a, b, c) of polynomials reduced modulo h(x). We assume h is a kernel polynomial, so alpha represents the restriction of an isogeny the kernel cut out by h.
    h: kernel polynomial for an ell isogeny for some odd prime ell
    """
    (a, b, c) = alpha
    (u, s, d) = phi
    max_deg = max(u.degree(),s.degree(),d.degree())
    a_h = [1, a]
    c_h = [1, c]
    while len(a_h) <= max_deg:
        a_h.append(a_h[-1] * a % h)
        c_h.append(c_h[-1] * c % h)
    ua = sum(u[i] * a_h[i] % h * c_h[u.degree()-i] % h for i in range(u.degree() + 1))
    sa = sum(s[i] * a_h[i] % h * c_h[s.degree()-i] % h for i in range(s.degree() + 1))
    da = sum(d[i] * a_h[i] % h * c_h[d.degree()-i] % h  for i in range(d.degree() + 1))
    sa = sa * b % h 
    da = da * c % h 
    return (ua, sa, da)

def proj_compose_and_reduce_chain(chain, h, phi=None):
    """Returns the coordinate functions of the isogeny phi represented by chain,
     a list of isogenies, modulo h. phi is an optional input which is
     precomposed with chain."""
    if phi:
        (a, b, c) = phi
    else:
        x = h.parent().gen()
        (a, b, c) = (x, 1, 1)
    for psi in chain:
        (a, b, c) = proj_isogeny_compose_mod(psi, (a, b, c), h)
    return (a, b, c)

def proj_trace_of_endo_mod(E, chain, deg, ell):
    """Computes the trace of endo mod ell. We require gcd(ell, deg) = 1. 

    endo: an endomomorphism of a supersingular elliptic curve E
    deg: the degree of endo
    ell: an odd prime
    
    This could be generalized to the case gcd(ell, deg) != 1: first compute any 
    kernel polynomial h defining ker phi. If the computation fails, then ker phi 
    is contained in ker endo. Compute a second h; if the computation succeeds, 
    then the output is correct. Otherwise, E[ell] is contained in ker endo, endo 
    factors through [ell], and its trace is 0 modulo ell. 
    """
    F.<x> = PolynomialRing(E.base_ring())
    A = E.a4()
    B = E.a6()                     
    h = get_kernel_poly(E, ell)
    f = x^3 + A*x + B
    endo_ell = proj_compose_and_reduce_chain(chain, h)
    endo_ell_2 = proj_compose_and_reduce_chain(chain, h, endo_ell)
    deg_ell = scale(deg % ell, (x, 1, 1), A, f, h)
    # compute characteristic equation endo_ell_2 + deg_ell = t * endo_ell
    LHS = proj_add(endo_ell_2, deg_ell, A, f, h)
    if iszero(LHS, h): return 0
    if isequal(LHS, endo_ell, h): return 1
    if isequal(LHS, negate(endo_ell, h), h): return -1
    RHS = endo_ell
    for t in range(2, ell - 1):
        RHS = proj_add(RHS, endo_ell, A, f, h)
        if isequal(LHS, RHS, h):
            return t

def trace_of_endo_mod_char(endo):
    """
    Computes the trace of the endomorphism modulo the characteristic p
    as the field trace of the scaling factor of the endomorphism
    """
    u = endo.scaling_factor()
    f = u.minpoly()
    assert f.degree() > 0
    if f.degree() == 1:
        f = f^2  # charpoly = minpoly^2
    t = (-f[1]).lift_centered()
    return t

def superSEA(endo):
    """Computes the trace of endo.
    
    endo: a separable endomorphism of a supersingular elliptic curve."""
    E = endo.domain()
    deg = endo.degree()
    proj_chain = projective_chain(endo)
    t, N = 0, 1
    p = E.base_field().characteristic()
    tp = trace_of_endo_mod_char(endo)
    t = crt(t, tp, N, p)
    N = lcm(N, p)
    ell = 2
    bound = 4 * sqrt(deg)
    while N <= bound:
        ell = next_prime(ell)
        if ell == p:
            ell = next_prime(ell)
        while deg % ell == 0:
            ell = next_prime(ell)
        tell = proj_trace_of_endo_mod(E, proj_chain, deg, ell)
        a = N * N.inverse_mod(ell)
        b = ell * ell.inverse_mod(N)
        N *= ell
        t = (a*tell + b*t) % N
    if t >= N/2: return t - N
    else: return t

