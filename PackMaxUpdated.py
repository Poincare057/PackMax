import numpy as np
import csv

##BEGIN PARAMETERS OF PROBLEM SIZE
prange = 37
irange = 16
##END PARAMETERS OF PROBLEM SIZE


##BEGIN RELEVANT ARRAYS
A = np.zeros([prange, 6])
T = np.zeros([irange, 5])
F = np.array([1, 4, 5, 2, 3, 4, 5, 6, 1, 1, 1, 3, 4, 3, 2, 2])        #default values, will modify for user input later
##END RELEVANT ARRAYS


##BEGIN INITIALIZE A, T FROM CSV FILES
#UNITS: LENGTH - DECIMETERS, MASS - KILOGRAMS
prod = open('PackMaxProducts.csv', 'r')
prod_reader = csv.reader(prod)
products = []
for line in prod_reader:
    products += [line]
for p in range(len(products) - 2):
    A[p][0] = int(round(float(products[p + 2][8])/100))
    A[p][1] = int(round(float(products[p + 2][7])/100))
    A[p][2] = int(round(float(products[p + 2][6])/100))
    A[p][3] = int(round(float(products[p + 2][5])))
    if p in range(26, 31):
        A[p][4] = int(products[p + 2][10])
    else:
        A[p][4] = int(products[p + 2][9])
    A[p][5] = 4 #default value, will modify for user input later
prod.close()

#trucks
truck = open('PackMaxTrucks.csv', 'r')
truck_reader = csv.reader(truck)
trucks = []
for line in truck_reader:
    trucks += [line]
for i in range(len(trucks) - 1):
    T[i][0] = int(round(3.048*float(trucks[i + 1][1])))
    T[i][1] = int(round(3.048*float(trucks[i + 1][2])))
    if i < 15:
        T[i][2] = int(round(3.048*float(trucks[i + 1][3])))
    else:
        T[i][2] = int(round(3.048*8))
    T[i][3] = int(round(1000*float(trucks[i + 1][4][:-3])))
    T[i][4] = 2 #default value, will modift for user input later
truck.close()
##END INITIALIZE A, T FROM CSV FILES


##BEGIN MORE RELEVANT ARRAYS
jrange = 0 
crange = prange
crange1 = 0 
for i in range(irange):
    if T[i][4] > jrange:
        jrange = T[i][4]
    crange += 2*T[i][4]
for p in range(prange):
    for i in range(irange):
        for j in range(int(T[i][4])):
            crange1 += 1
jrange, crange = int(jrange), int(crange) + crange1

M = np.zeros([prange, irange, jrange])
gradM = np.zeros([prange, irange, jrange])
C = np.zeros([crange, prange, irange, jrange])
B = np.zeros([crange])
##END MORE RELEVANT ARRAYS


##BEGIN AUXILIARY FUNCTIONS
def r(p, i, j):
    return(1 - 1000*F[i]*A[p][3]/T[i][3])

def dot_compare(C, M, B):
    res = True
    for c in range(crange):
        dot = 0
        for p in range(prange):
            for i in range(irange):
                for j in range(int(T[i][4])):
                    dot += C[c][p][i][j]*M[p][i][j]
        if dot > B[c]:
            res = False
            break
    return(res)

def displayM(M):
    for p in range(prange):
        for i in range(irange):
            for j in range(int(T[i][4])):
                print(p, i, j, M[p][i][j])
        print()
        
##END AUXILIARY FUNCTIONS


##BEGIN CONSTRAINTS
#CONSERVATION OF NUMBER
temp_rank = 0
for p in range(prange):
    for i in range(irange):
        for j in range(int(T[i][4])):
            C[temp_rank][p][i][j] = 1
    B[temp_rank] = A[p][5]
    temp_rank += 1

#NOT EXCEEDING TONNAGE
for i in range(irange):
    for j in range(int(T[i][4])):
        for p in range(prange):
            C[temp_rank][p][i][j] = A[p][3]*A[p][4]
        B[temp_rank] = T[i][3]
        temp_rank += 1

maxd, maxw, maxh = 0, 0, 0 
for p in range(prange):
    if A[p][0] > maxd:
        maxd = A[p][0]
    if A[p][1] > maxw:
        maxw = A[p][1]
    if A[p][2] > maxh:
        maxh = A[p][2]

#PACKABILITY
for i in range(irange):
    for j in range(int(T[i][4])):
        for p in range(prange):
            C[temp_rank][p][i][j] = maxd*maxw*maxh
        B[temp_rank] = T[i][0]*T[i][1]*T[i][2]

#NONNEGATIVITY
for i in range(irange):
    for j in range(int(T[i][4])):
        for p in range(prange):
            C[temp_rank][p][i][j] = -1
            B[temp_rank] = 0
            temp_rank += 1
##END CONSTRAINTS


##BEGIN SOLVER
for i in range(irange):
    for j in range(int(T[i][4])):
        for p in range(prange):
            gradM[p][i][j] = r(p, i, j)

gamma = 0.1
t = 2
while gamma > 0.0001:
    dcom = dot_compare(C, M, B)
    if dcom:
        M = np.round(M + gamma*gradM)
        #M = M + gamma*gradM
    else:
        M = np.round(M - gamma*gradM)
        #M = M - gamma*gradM
    t += 1
    #gamma = gamma*(t/(t + 1))
    gamma = 1/t**2
    print(dcom)         #debug

##END SOLVER


##BEGIN PACKING EACH TRUCK


##END PACKING EACH TRUCK
