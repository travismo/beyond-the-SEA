# beyond-the-SEA
Sage code to compute the trace of a separable endomorphism of a supersingular elliptic curve. The file cycles.sage also includes code for computing endomorphisms of supersingular elliptic curves via the algorithm of https://arxiv.org/abs/2004.11495. The implementation in superSEA.sage is loosely based on Andrew Sutherland's implementation of Schoof's algorithm found here: https://cocalc.com/AndrewVSutherland/18.783EllipticCurves2023/SchoofsAlgorithm.

We assume the input to our algorithm is an endomorphism of a supersingular elliptic curve E over a finite field Fp2 of p^2 elements whose j-invariant is not 0 or 1728. We exploit the fact that, given our assumptions on E, every ell isogeny of E is defined over Fp2 -- every prime is an "Elkies prime" for E! Thus, for the first sufficiently many odd primes ell, we compute the trace of the endomorphism modulo ell by restricting the endomorphism to the kernel cut out by the kernel polynomial of some ell-isogeny. By working with projective coordinates, we avoid inversions in polynomial quotient rings. Finally, we use the fact that we can compute the trace of the endomorphism modulo p using the "scaling factor" of the endomorphism. 

Requires Sage 10.5.

Sample usage:
```
load('cycles.sage')
load('isogenies.sage')
load('superSEA.sage')
p = next_prime(2**16)
# find one supersingular j-invariant
j = find_supersingular_j(p)
# find a pseudorandom supersingular j-invariant as the final vertex in a random walk in the 2-isogeny graph
j = random_walk(j,2)[-1]
# compute a cycle in the 2-isogeny graph based at j
C = CGL_cycle(j,2)
# convert the cycle into an endomorphism of a curve with j-invariant j
endo = path_to_isogeny(C, 2)
# compute the trace of the endomorphism
t = superSEA(endo)
```

Authors: Travis Morrison, Lorenz Panny, Jana Sotakova, Michael Wills. 


