#!/usr/bin/env python3

from pprint      import pprint
from pyeda.inter import exprvar, exprvars, expr2bdd, expr, bddvar
from itertools   import accumulate, product

def inMain(f):
    if __name__ == '__main__':
        f()

def printHeader(string):
    end = '-'*20
    print('{} {:^20} {}'.format(end,string,end))

def showDictionary(d):
    return '\n'.join(d.map(lambda x,y: '{x}: {y}'.format(x,y)))

def getBin(i, length):
    return [ int(i) for i in str(bin(i)[2:]).zfill(length) ]

# step 1

printHeader("Step 1")

G = [ (i,getBin(i, 5)) for i in range(32) ]

print("G: [(number, [bits])]")
pprint(G)

# step 2

printHeader("Step 2")

G_Edges = [ iv + jv for (i, iv) in G for (j, jv) in G
            if (i+3)%32==j%32 or (i+8)%32==j%32]

print("G_Edges: [ 5 bits of i + 5 bits of j ]")
pprint(G_Edges)

primes = [3,5,7,11,13,17,19,23,29,31]
evens  = [i for i in range(32) if i % 2 == 0]
print("primes:")
pprint([ getBin(i, 5) for i in primes])
print("evens:")
pprint([ getBin(i, 5) for i in evens])

# step 3

printHeader("Step 3")

def parseVariable(name, value):
    return '~' + name if value == 0 else name

names = [ i+j for (i,j) in product('xy', '12345') ]
def transformEdge(edge, names):
    values = [ parseVariable(v,i) for i,v in zip(edge, names) ]
    bdd = ' & '.join(values)
    return bdd

bdds = map(lambda x: transformEdge(x, names), G_Edges)

print("made bdds")

# step 4

printHeader("Step 4")

def join_disjunct(a):
    return '(' + ') | ('.join(a) + ')'

bdds_together = join_disjunct(bdds)

print("bdds_together: long string bdd for entire equation")
pprint(bdds_together)

# step 4'

printHeader("Step 4'")

x_names = names[0: 5]

p_bdds = join_disjunct(map(lambda x: transformEdge(x, x_names), map(lambda x: getBin(x, 5), primes)))
p = expr2bdd(expr(p_bdds))

y_names = names[5:-1]
e_bdds = join_disjunct(map(lambda y: transformEdge(y, y_names), map(lambda y: getBin(y, 5), evens)))
e = expr2bdd(expr(e_bdds))

print("p_bdds: like bdds_together")
pprint(p_bdds)
print("e_bdds: like bdds_together")
pprint(e_bdds)

# step 5

printHeader("Step 5")

x1, x2, x3, x4, x5, y1, y2, y3, y4, y5 = map(bddvar, names)
r = expr2bdd(expr(bdds_together))
xx1, xx2, xx3, xx4, xx5 = [ bddvar('xx' + str(i + 1)) for i in range(5) ]
yy1, yy2, yy3, yy4, yy5 = [ bddvar('yy' + str(i + 1)) for i in range(5) ]
def rr(r, xx1, xx2, xx3, xx4, xx5, yy1, yy2, yy3, yy4, yy5):
    return r.compose({x1: xx1, x2: xx2, x3: xx3, x4: xx4, x5: xx5,
                      y1: yy1, y2: yy2, y3: yy3, y4: yy4, y5: yy5})

print("Made expression Tree and made variables")

# step 5'

printHeader("Step 5'")
pp = rr(p, xx1, xx2, xx3, xx4, xx5, xx1, xx2, xx3, xx4, xx5)
ee = rr(e, yy1, yy2, yy3, yy4, yy5, yy1, yy2, yy3, yy4, yy5)

print ("Made expression trees for primes and evens")

# step 6

printHeader("Step 6")
zz1, zz2, zz3, zz4, zz5 = [ bddvar('zz' + str(i + 1)) for i in range(5) ]

rrxz = rr(r, xx1,xx2,xx3,xx4,xx5,zz1,zz2,zz3,zz4,zz5)
rryz = rr(r, yy1,yy2,yy3,yy4,yy5,zz1,zz2,zz3,zz4,zz5)
rr2 = (rrxz & rryz).smoothing([zz1, zz2, zz3, zz4, zz5])

print("Made RR2")

# step 7
printHeader("Step 7")

hh  = rr2
hhp = rr2 & False
while not hhp.equivalent(hh):
    hhp  = hh
    hh   = hh | (hh & rr2).smoothing([zz1, zz2, zz3, zz4, zz5])
    rr2p = hh

print("Made rr2p")

# step 8

printHeader("Step 8")
jj = (rr2p & ee).smoothing([yy1, yy2, yy3, yy4, yy5])
qq = (~(~(jj | ~pp))).smoothing([xx1, xx2, xx3, xx4, xx5])

print("QQ equivalent to bdd true? ", bool(qq))
