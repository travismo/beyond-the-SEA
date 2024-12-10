def normalized_model(E, j_tilde, ell):
    """Computes an elliptic curve E~ with the given j-invariant such that there
    exists an ell isogeny from E to E~.

    If E has an ell-isogeny to an elliptic curve E~ with j-invariant j_tilde,
    then this function will return the coefficients A_tilde and B_tilde such
    that E~ is defined by y^2 = x^3 + A_tilde*x + B_tilde. This could be extended to also output the sum of the x-coordinates of the kernel of the ell-isogeny.

    Inputs:
        E: An elliptic curve defined over the finite field F of characteristic p. Must have j(E) not in {0, 1728}
        j_tilde: An element of F. Must have that (j(E), j_tilde) is a single root of the ell'th modular polynomial
        ell: A positive integer with ell + 2 < p if p > 0. 
    Outputs:
        A_tilde, B_tilde: Elements of F
        """
    # Set-up
    j = E.j_invariant()
    A = E.a_invariants()[3]
    B = E.a_invariants()[4]
    FXY.<X,Y> = PolynomialRing(E.base_field())

    if j == 0:
        raise NotImplementedError('j = 0')
    if j == 1728:
        raise NotImplementedError('j = 1728')

    # Compute all deriviatives and partial derivatives of the ell'th modular polynomial at (j, j_tilde)
    # PHI = ClassicalModularPolynomialDatabase()
    Phi_ell = classical_modular_polynomial(ell)
    Phi_ell = FXY(Phi_ell)
    Phi_X = Phi_ell.derivative(X)
    Phi_Y = Phi_ell.derivative(Y)
    # Phi_XX = Phi_X.derivative(X)(j,j_tilde)
    # Phi_YY = Phi_Y.derivative(Y)(j,j_tilde)
    # Phi_XY = Phi_X.derivative(Y)(j,j_tilde)
    Phi_X = Phi_X(j,j_tilde)
    Phi_Y = Phi_Y(j,j_tilde)

    # Recover the coefficients of E~
    m = 18*B/A
    j_prime = m*j
    # k = j_prime / (1728-j)
    j_tilde_prime = -j_prime * Phi_X / (ell * Phi_Y)
    m_tilde = j_tilde_prime / j_tilde
    k_tilde = j_tilde_prime / (1728 - j_tilde)
    A_tilde = ell^4 * m_tilde * k_tilde / 48
    B_tilde = ell^6 * m_tilde^2 * k_tilde / 864

    # Recover the sum of the abscissas of the kernel
    # r = -(j_prime^2*Phi_XX
    #       + 2*ell*j_prime*j_tilde_prime*Phi_XY
    #       + ell^2*j_tilde_prime^2*Phi_YY) / (j_prime*Phi_X)
    # p1 = ell*(r/2 + (k-ell*k_tilde)/4 + (ell*m_tilde - m)/3)

    return A_tilde, B_tilde #, p1


def get_kernel_poly(E, ell):
    """Computes a kernel polynomial for an ell-isogeny of E.
    
    E: a supersingular elliptic curve
    ell: an odd prime"""
    j = E.j_invariant()
    F = classical_modular_polynomial(ell, j)
    x = F.parent().gen()
    while F.degree() > 0:
        j1 = F.any_root()
        if F.derivative()(j1) != 0 and j1 not in [F(0), F(1728)]:
            A_tilde, B_tilde = normalized_model(E, j1, ell)
            Etilde = EllipticCurve([A_tilde, B_tilde])
            h = compute_isogeny_bmss(E, Etilde, ell)
            return h
        F = F // (x - j1)^F.valuation(x - j1)
    # all roots of Phi(j, x) have multiplicity
    return E.isogenies_prime_degree(ell)[0].kernel_polynomial()



    



