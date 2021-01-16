from sympy import *

# Originally:
#
# rho_ij,kl := basis_L * rho * basis_R
#
# In each index pair, the first index indexes the basis of the first system,
# while the second index indexes the basis of the second system.
#
# In partial transpose:  ij,kl -> il,kj
#
def partial_transpose(rho):
    from sympy.physics.quantum import TensorProduct as kron

    # Standard orthonormal basis.
    ex = Matrix( [S("1"), S("0")] )  # Matrix creates column vectors by default
    ey = Matrix( [S("0"), S("1")] )

    # Bases for first and second systems.
    #
    e1 = ( ex, ey )
    e2 = ( ex, ey )

    # map index pair to combined index
    # (same ordering of basis elements as in Example II.12 on p. 7)
    p = { (0,0) : 0, (0,1) : 1, (1,0) : 2, (1,1): 3 }

    # ij,kl -> il,kj
    n = len(e1)  # assumed same as len(e2)
    rho_pt = zeros(2*n)
    for i in range(n):
        for j in range(n):
            basis_L = kron( e1[i], e2[j] ).adjoint()  # adjoint() = dagger
            for k in range(n):
                for l in range(n):
                    basis_R = kron( e1[k], e2[l] )
                    rho_pt[ p[(i,l)], p[(k,j)] ] = basis_L * rho * basis_R

    return rho_pt

# Same as above, but hardcoded for the standard basis, which allows us
# to implement this by simply swapping indices.
#
def partial_transpose_simple(rho):
    n = 2
    p = { (0,0) : 0, (0,1) : 1, (1,0) : 2, (1,1): 3 }

    rho_pt = zeros(2*n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    rho_pt[ p[(i,l)], p[(k,j)] ] = rho[ p[(i,j)], p[(k,l)] ]

    return rho_pt


def main():
    # build rho and compute its eigenvalues
    #
    p = S("p")
    a = p/2
    b = (1-p)/2
    rho = Matrix( [[a, 0, 0, a], [0, b, b, 0], [0, b, b, 0], [a, 0, 0, a]] )

    a = p/2
    b = (1+p)/4
    c = (1-p)/4
    rho1 = Matrix( [[b, 0, 0, a], [0, c, 0, 0], [0, 0, c, 0], [a, 0, 0, b]] )

    for r in [rho, rho1]:
        print "=" * 80
        pretty_print(r)
        print r.eigenvals()

        rho_pt = partial_transpose(r)
        pretty_print(rho_pt)
        #pretty_print( partial_transpose_simple(r) )  # same as rho_pt
        print rho_pt.eigenvals()
        pretty_print(rho_pt.eigenvects())


if __name__ == '__main__':
    main()

