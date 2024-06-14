#https://arxiv.org/abs/2309.13741
BasisOrd = {} #given basis element, give the order
BasisEle = [] #contains basis elements
tmp = []
nn = 0 #dimension of input matrix
kk = 0
TT = []
XX = []
def ord(X): #Orders the basis
    r = 0
    for i in range(0,len(X)):
        if X[i]!=0:
            r+=1
    return r
def go(b,c):
    if b==kk:
        #print("wtf ",c)
        TT.append(c)
        return
    for aa in range(0,nn):
        if XX[aa]!=0:
            cc = [c[d] for d in range(0,len(c))]
            cc.append(aa)
            XX[aa]-=1
            go(b+1,cc)
            XX[aa]+=1
def CreateOrders(X):
    global TT,XX
    TT = []
    XX = X
    go(0,[])
    return TT
def hagale(a,b,s): #this creates all relevant tuples (index of basis)
    global BasisEle,BasisOrd,tmp
    if b==kk+1:
        yes = 0
        for z in range(0,nn):
            if tmp[z]==kk:
                yes = 1
            if yes==1:
                break
        if yes==1:
            return
        BasisEle.append(tmp.copy())
        BasisOrd[s[1:]]=len(BasisEle)
        #print(tmp,s[1:],normBasis(tmp),BasisEle)
        return
    if a==nn:
        tmp[a-1]+=1
        ss = s+","+str(a)
        hagale(a,b+1,ss)
        tmp[a-1]-=1
        return
    tmp[a-1]+=1
    ss = s+","+str(a)
    hagale(a,b+1,ss)
    tmp[a-1]-=1
    hagale(a+1,b,s)
def normBasis(X): #this is the norm of orthogonal basis, (1/this) normalizes norm
        s = factorial(kk)
        for i in range(0,nn):
            s/=factorial(X[i])
        return sqrt(s)
def CreateBasis(n,k): #we create the basis here
    global BasisOrd,BasisEle, tmp,nn,kk
    nn = n
    kk = k
    BasisOrd = {}
    BasisEle = []
    for i in range(0,n):
        ttt = [0 for ii in range(0,nn)]
        ttt[i]=k
        BasisEle.append(ttt)
        st = ""
        for jj in range(0,k):
            st+=","+str(i+1)
        BasisOrd[st[1:]]=i+1
    tmp = [0 for i in range(0,n)]
    hagale(1,1,"")
    BasisOrd=sorted(BasisOrd,key=ord)
def Coeff(A,i,j): #Computes coefficient (i,j) of big matrix USING MISTERIOUS COMBINATORICS

    co2 = 0
    x = BasisEle[i]
    xx = [0 for a in range(0,kk)]
    ll = 0
    #gives ordered sequence from indexed sequence on i
    for a in range(0,nn):
        for b in range(0,x[a]):
            xx[ll]=a
            ll+=1
    y = BasisEle[j]
    co =  normBasis(x)/normBasis(y)
    #takes care of order on i
    #for a in range(0,nn):
    #    co*=factorial(x[a])
    #iterate over orderings on j (all orderings?) multiply A[ordered index on x][index on ordering of y]
    L = CreateOrders(y)
    #print("Creating orders ",i,j,y,L)
    for a in range(0,len(L)):
        co3 = 1
        for b in range(0,kk):
            co3*=A[L[a][b]][xx[b]]
        co2+=co3
    return co*co2
def Coeff2(A,i,j): #Computes coefficient (i,j) of big matrix USING THE INNER PRODUCT
    co2 = 0
    x = BasisEle[i]
    xx = [0 for a in range(0,kk)]
    ll = 0
    #gives ordered sequence from indexed sequence on i
    for a in range(0,nn):
        for b in range(0,x[a]):
            xx[ll]=a
            ll+=1
    y = BasisEle[j]
    co =  1/(normBasis(y)*normBasis(x))
    #takes care of order on i
    #for a in range(0,nn):
    #    co*=factorial(x[a])
    #iterate over orderings on j (all orderings?) multiply A[ordered index on x][index on ordering of y]
    L1 = CreateOrders(y)
    L2 = CreateOrders(x)
    #print("Creating orders ",i,j,y,L)
    for a in range(0,len(L1)):
        for c in range(0,len(L2)):
            co3 = 1
            for b in range(0,kk):
                co3*=A[L1[a][b]][L2[c][b]]
            co2+=co3
    return co*co2
#given permutation p, generates matrix represent it
def matToPerm(A):
    AA = []
    for i in range(0,len(A)):
        for j in range(0,len(A)):
            if A[i][j]==1:
                AA.append(j)
                break
    return AA
def permToMat(p):
    A = []
    for i in range(0,len(p)):
        A.append([])
        for j in range(0,len(p)):
            if p[i]==j:
                A[i].append(1)
            else:
                A[i].append(0)
    return A
#generatesLineGraph
def genLineGraph(n):
    AA = []
    for i in range(0,n):
        AA.append([])
        for j in range(0,n):
            if j==i+1 or j==i-1:
                AA[i].append(1)
            else:
                AA[i].append(0)
    return AA
def Symk(A): #computes A**Symk
    AA = [[Coeff2(A,j,i) for j in range(0,len(BasisEle))] for i in range(0,len(BasisEle))]
    return AA
def CreateGraph(A,di):
    G = Graph(loops=true)
    if di:
        G = DiGraph(loops=true)
    for i in range(0,len(A)):
        G.add_vertex(i)
    for i in range(0,len(A)):
        for j in range(0,len(A)):
            if A[i][j]!=0:
                G.add_edge((i,j),label=str(A[i][j]))
                G.set_edge_label(i,j,str(A[i][j]))
    return G
#ShowGraph muestra el grafo, tienes que meter en A la matriz de adyacencia (en lista de lista, como A = [[1,0],[0,1]]),
#en k el exponente del Symk, en s el tamagno de la imagen(numero entero), en ed True si quieres ver el peso de las aristas o False si no,
# en di True si el grafo es dirigido o False si no.
def ShowGraph(A,k,s,ed,di):
    CreateBasis(len(A),k)
    AA = Symk(A)
    #print(matrix(AA))
    CreateGraph(AA,di).plot(edge_labels=ed).show(figsize = [s,s])
# Dada una dimension y k, genera una matriz con variables aij y computa la matriz A^Sym_k
def printGenMatrix(n,k):
    A = [[0 for j in range(0,n)] for i in range(0,n)]
    for i in range(0,n):
        for j in range(0,n):
            A[i][j]=var("a"+str(i+1)+str(j+1))
    CreateBasis(n,k)
    AA = Symk(A)
    print(matrix(AA))
    return AA
#Imprime la matriz A^Sym_k.
def printSymk(A,k):
    CreateBasis(len(A),k)
    AA = Symk(A)
    print(matrix(AA))
    return AA
print("Put adjacency matrix and k and computes matix ^Sym_k after pressing enter")
@interact
def _(A = input_box('[[1,1],[1,0]]',type= list ,label='Matrix'),k=input_box('3',type= int ,label='k'),s=input_box('5',type= int ,label='size to draw'),ed = checkbox(False, "Show edge weight?"),di= checkbox(False, "Directed?")):
    ShowGraph(A,k,s,ed,di)
