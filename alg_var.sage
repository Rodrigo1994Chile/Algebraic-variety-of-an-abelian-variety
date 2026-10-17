def decomposing(f):
    """
    INPUT
    =====
    f               :
        Polynomial in variables [x1,...,xn] with coefficients in a field extension K/Q such that K=<a1,...,ar> as Q-vector space

    OUTPUT
    ======
    Fi=[F1,...,Fr]  :
        List of polynomials in variables [x1,...,xn] with coefficients in Q such that f= F1*a1+...+Fr*ar
    """
    R=f.parent()
    K=R.base_ring()
    n=K.degree()
    A=R.gens()
    m=len(A)
    Fi=[0 for i in range(n)]
    for r in range(m):
        mon=f[A[r]]
        for i in range(n):
            aux=Fi[i]
            aux+=mon[i]*A[r]
            Fi[i]=aux
    return Fi

def alg_var(d,Z,n):
    """
    INPUT
    =====
        d    :
            List of positive integers d1,...,dg such that d1|...|dg
        Z       :
            Matrix of gxg in the Siegel upper half space
        n       :
            Positive integer between 1 and g-1
    OUTPUT
    =====

        polynomial_equations :
            List of lists whose first element is a polynomial in variables (a_{i,j}) with i<j<g and the second element is the right hand side of the equations
    """
    # ============================================================
    # Constants throughout the routine
    # ============================================================
    
    g=len(d)
    R=Z[0,0].parent()
    polynomial_equations=[]
    
    # ============================================================
    # Define the variables of the algebraic variety
    # ============================================================
    
    nam=[]
    for i in range(1,2*g):
        for j in range(i+1,2*g+1):
            nam.append("a"+str(i)+"c"+str(j))

    # ============================================================
    # Operations work with polynomials in the variables
    #    previously define
    # ============================================================
    
    R=PolynomialRing(R,len(nam),nam)
    a=R.gens()

    # ============================================================
    # Set the diagonal matrix of the polarization type and the 
    #     period matrix Pi=[D Z]
    # ============================================================
    D=zero_matrix(R,g,g)
    for i in range(g):
        D[i,i]=d[i]
    Pi=block_matrix([[D,Z]])

    # ============================================================
    # Set the exterior algebra of (real) 1-forms
    # ============================================================
    E=ExteriorAlgebra(R,names=tuple( ['dx'+str(i) for i in range(1, 2*g+1)]))
    dx=E.gens()

    # ============================================================
    # Define the complex 1-forms with the period matrix 
    #      according to theorem
    # ============================================================
    
    dz=[]
    for i in range(g):
        aux=0
        for j in range(2*g):
            aux+=Pi[i,j]*dx[j]
            aux
        dz.append(aux)
        
    # ============================================================
    # Define a generic 1-form in the manifold
    # ============================================================
    
    omega=0
    cont=0
    for i in range(1,2*g):
        for j in range(i+1,2*g+1):
            omega+=a[cont]*dx[i-1]*dx[j-1]
            cont+=1

    # ============================================================
    # Define the 1-form corresponding to the canonical polarization
    # ============================================================
    
    w0=0
    for k in range(g):
        w0+=d[k]*dx[k]*dx[k+g]
    w0=-w0
    w0

    # ============================================================
    # Get the (rational) linear equations defining Y_Pi
    # ============================================================
    
    aux=omega
    for i in range(g):
        aux=aux*dz[i]
    for mon in aux.terms():
        for rational_eq in decomposing(mon.coefficients()[0]):
            polynomial_equations.append([rational_eq,0])
            
    # ============================================================
    # Get the non-linear equations defining V_Pi
    # ============================================================
    
    gamma=(-1)^g
    for i in range(g):
        gamma=gamma*dx[i]*dx[i+g]
    X=1
    for i in range(g):
        X*=d[i]
    for r in range(1,g+1):
        polynomial_equations.append([])
        polynomial_equations[-1].append((omega^r*w0^(g-r)).terms()[0].coefficients()[0])
        if n<r:
            polynomial_equations[-1].append(0)
        else:
            aux=1
            for k in range(n+1,g+1):
                aux*=k-r
            polynomial_equations[-1].append(X*factorial(n)*aux)

    # ============================================================
    # Return the list of polynomial equations
    # ============================================================
    
    return polynomial_equations
