import numpy as np
import copy
import math
import random
import pylab as pl


#DEB Dataset
class Individual():
    def __init__(self,x):
        self.x=x
        self.NumX=len(x)
        self.crowd=0.0
        self.d=0   #dominate
        self.bd=0  #by dominated
        f1=float(x[0]/1000)
        h=float(1+x[1]/100)
        f2=h*(1-(f1/h)**2-(f1/h)*np.sin(8*np.pi*f1))
        self.f=[f1,f2]    #multiobjective function
#cal x dominated y or not
def Dominate(x,y,min=True):
    if min:

        for i in range(len(x.f)):
            if x.f[i]>y.f[i]:
                return False


        return True
    else:
        for i in range(len(x.f)):
            if x.f[i]<y.f[i]:
                return False
        return True
def UpdateDomination(B,Nd):
    #B:pre updated population,Nd :maximum number of domination population
    #output Nd dominated solutions
    n=len(B)
    DT=[]
    for i in range(n):
        for j in range(i+1,n):
            if Dominate(B[i],B[j]):
                B[i].d+=1
                B[j].bd+=1
    for i in range(n):
        if B[i].bd==0:
            DT.append(B[i])
    if len(DT)<=Nd:
        return DT
    else:
        DT.sort(key=lambda p:p.crowd,reverse=True)
        return DT[:Nd]
def ActiveSelection(D,na):
    if len(D)<=na:
        return D
    else:
        D.sort(key=lambda p:p.crowd,reverse=True)
        return D[:na]

#cal every individual's crowd distance
def CrowdDistance(B):
    n=len(B)
    for i in range(len(B[0].f)):
        B.sort(key=lambda p:p.f[i])   #sorted by values of f1 and f2
        B[0].crowd=float('inf')
        B[n-1].crowd=float('inf')
        h=B[n-1].f[i]-B[0].f[i]
        for j in range(1,n-1):
            B[j].crowd+=((B[j+1].f[i]-B[j-1].f[i])/h)
def Clone(A,Nc):
    #proportional clone operaton
    num=0.0
    c=[]
    A.sort(key=lambda p:p.crowd,reverse=True)
    for i in range(len(A)):
        if A[i].crowd!=float('inf'):
            k=i
            break
    for j in range(k):
        A[j].crowd=A[k].crowd*2


    for i in range(len(A)):
        num+=A[i].crowd
    for i in range(len(A)):
        l=int((A[i].crowd*Nc)/num)
        for j in range(l):
            c.append(A[i])
    return c
def CrossOver(C,A):
    T=[]
    for i in range(len(C)):
        r=random.randint(0,len(A)-1)
        a=A[r]
        Xa=C[i].x
        Xb=a.x
        xtemp = []
        for j in range(len(a.x)):
            str1,str2=integerToString(int(Xa[j]),int(Xb[j]))
            R=random.randint(0,len(str1)-1)
            temp1=str1[:R]+str2[R:]
            temp2=str2[:R]+str1[R:]
            xtemp.append(stringToInteger(temp1))
        T.append(Individual(xtemp))
    return T
def Mutate(T,pm):
    R=[]
    for i in range(len(T)):
        if random.random()<pm:
            str1,str2=integerToString(int(T[i].x[0]),int(T[i].x[1]))
            strc=str1+str2
            r=random.randint(0,len(strc)-2)
            if strc[r]=='0':
                strctempt=strc[:r]+'1'+strc[r+1:]
            else:
                strctempt = strc[:r] + '0' + strc[r + 1:]
            num1=stringToInteger(strctempt[:len(str1)])
            num2=stringToInteger(strctempt[len(str1):])
            while num1>=1000:
                num1=num1/10
            while num2>=1000:
                num2=num2/10
            R.append(Individual([num1,num2]))
        else:
            R.append(T[i])
    return R




def integerToString(numero,numero2):
    '''
        int numbers transformed to Binary as String

    '''
    cadenas = []
    cadenas.append(bin(numero))
    cadenas.append(bin(numero2))
    #Limpiar las cadenas
    cadenas[0] = cadenas[0].replace("0b","")
    cadenas[1] = cadenas[1].replace("0b","")
    if (len(cadenas[0]) > len(cadenas[1])):
        cadenas[1] = ("0"*(len(cadenas[0])-len(cadenas[1]))) + cadenas[1]
    else:
        cadenas[0] = ("0"*(len(cadenas[1])-len(cadenas[0]))) + cadenas[0]
    return cadenas[0], cadenas[1]

def stringToInteger(cadena):
    '''
       transformed Binary to Integer
                '''
    decimal = 0
    for i,v in enumerate(cadena):
        if(v == '1'):
            decimal = decimal + math.pow(2,len(cadena)-1-i)
    return decimal

def main(Gmax,nd,na,nc):
    #Gmax:maximum number of generations
    #nd:maximum size of dominant population)
    #na:maximum size of active population)
    #nc:size of clone population)
    B=[]
    # initial the population
    for i in range(nd):
        B.append(Individual([random.random()*1000,random.random()*1000]))
    for i in range(Gmax):
        print i
        CrowdDistance(B)
        D=UpdateDomination(B,nd)
        A=ActiveSelection(D,na)
        C=Clone(A,nc)
        CT=CrossOver(C,A)
        Ct=Mutate(CT,1.0/nc)
        B=Ct+D
    D = UpdateDomination(B, nd)
    x = []
    y = []
    for i in range(len(D)):
        x.append(D[i].f[0])
        y.append(D[i].f[1])
    pl.plot(x, y, '*')
    pl.xlabel('f1')
    pl.ylabel('f2')
    pl.show()
if __name__ == '__main__':
    main(200,200,50,200)




























