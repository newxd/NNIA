from numpy import *
import copy
import math
import random
import pylab as pl


#ZDT3 Dataset
class Individual():
    def __init__(self,x):
        self.x=x
        self.n=len(x)
        self.crowd=0.0

        self.d=0   #dominate
        self.bd=0  #by dominated as the rank
        f1=x[0]/1000
        num=0.0
        for i in range(1,len(x)):
            num+=(x[i]/1000)
        g=1.0+9.0*num/29
        #f2=g*(1-sqrt(f1/g)-f1/g*sin(10*pi*f1))
        f2=g*(1-(f1/g)**2)
        self.f=[f1,f2]    #multiobjective function
def Initial(Nd):
    #initial the population,number is Nd
    D=[]
    for i in range(Nd):
        x=[]
        for j in range(30):
            x.append(random.random()*1000)
        D.append(Individual(x))
    return D

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

def UpdateDomination(B,Nd,Na):
    #B:pre updated population
    # Nd :maximum number of domination population,Na<Nd
    #output  dominated solutions
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
    if len(DT)<Na:
        B.sort(key=lambda p:p.bd)
        DT=B[:Na]
        DSCD(DT,0.5)
        return DT

    elif len(DT)<=Nd:
        DSCD(DT,0.5)
        return DT
    else:
        DSCD(DT,0.5)
        DT.sort(key=lambda p:p.crowd,reverse=True)
        return DT[:Nd]

def ActiveSelection(D,na):
    if len(D)<=na:
        return D
    else:
        D.sort(key=lambda p:p.crowd,reverse=True)
        return D[:na]

def DSCD(B,alpha):
    #calculate every individual's Double-Sphere crowding distance
    n=len(B)
    B.sort(key=lambda p:p.f[0])
    B[0].crowd=float('inf')
    B[n-1].crowd=float('inf')
    for i in range(1,n-1):
        p1=NearestPoint(B[:i],B[i])
        p2=NearestPoint(B[i+1:],B[i])
        Ds1=(EuclideanDistance(p1,p2)/2)**2
        Ds2=((EuclideanDistance(B[i],p1)+EuclideanDistance(B[i],p2))/2)**2-Ds1
        B[i].crowd=alpha*Ds1+(1-alpha)*Ds2


def EuclideanDistance(a,b):
    #cal Euclidean distance between a and b
    num=0.0
    for i in range(len(a.f)):
        num+=(a.f[i]-b.f[i])**2
    return sqrt(num)


def NearestPoint(B,b):
    #find b's nearest individual in muliti_function space
    k=float('inf')
    for i in range(len(B)):
        num=0.0
        for j in range(len(b.f)):
            num+=(B[i].f[j]-b.f[j])**2
        if num<k:
            k=num
            best=B[i]
    return best

def RaClone(D,Nc):
    #RA Cloning operation,D:current population,Nc:clone size
    #output cloned population
    D.sort(key=lambda p:p.bd)
    rnk = []         #population's rank list
    temp = []
    r =D[0].bd
    temp.append(D[0])
    for i in range(1,len(D)):
        if D[i].bd!=r:
            rnk.append(temp)
            r=D[i].bd
            temp=[]
            temp.append(D[i])
        else:
            temp.append(D[i])
    rnk.append(temp)
    if len(rnk)<2:
        return Clone(D,Nc)
    Srank=[]                         #size of individuals in each rank
    for i in range(len(rnk)):
        Srank.append(float(Nc)/len(D)*len(rnk[i]))

    r=float(len(rnk[0]))/len(D)
    cloneP=[]
    if r<=0.33:     #early model
        a=float(len(rnk[0]))/(len(rnk[0])+len(rnk[1]))
        Srank[0]+=(1-a)*Srank[1]
        Srank[1]=a*Srank[1]
        for i in range(len(rnk)):
            cloneP.extend(Clone(rnk[i],int(Srank[i])))
        return cloneP
    elif r<=0.66:   #medium model
        num=0.0
        for i in range(1,len(Srank)):
            num+=Srank[i]
            Srank[i]=Srank[i]*r
        Srank[0]+=(1-r)*num
        for i in range(len(Srank)):
            cloneP.extend(Clone(rnk[i],int(Srank[i])))
        return cloneP
    else:          #late model
        for i in range(len(Srank)):
            cloneP.extend(Clone(rnk[i],int(Srank[i])))
        return cloneP

def Clone(B,Nc):
    #proportional clone operaton
    num=0.0
    c=[]
    k=-1
    A=sorted(B,key=lambda p:p.crowd,reverse=True)
    for i in range(len(A)):
        if A[i].crowd!=float('inf'):
            k=A[i].crowd
            break

    for i in range(len(A)):
        if A[i].crowd ==float('inf'):
            A[i].crowd = k * 2
        num += A[i].crowd
    if k==-1 or num==0:
        n=Nc/len(A)
        for i in range(len(A)):
            for j in range(int(n)):
                c.append(A[i])
    else:

        for i in range(len(A)):
            l=(A[i].crowd/num)*Nc
            for j in range(int(l)):
                c.append(A[i])
    return c

def CrossOver(C,A,pc):
    T=[]
    for i in range(len(C)):
        if random.random()<pc:
            r=random.randint(0,len(A)-1)
            a=A[r]
            Xa=C[i].x
            Xb=a.x
            xtemp = []
            for j in range(len(a.x)):
                str1,str2=integerToString(int(Xa[j]),int(Xb[j]))
                R=random.randint(0,len(str1)-1)
                #R=int(len(str1)/2)
                temp1=str1[:R]+str2[R:]
                temp2=str2[:R]+str1[R:]
                xtemp.append(stringToInteger(temp1))
            T.append(Individual(xtemp))
        else:
            T.append(C[i])
    return T

'''def Mutate(T,pm):
    R=[]
    n=len(T[0].x)
    for i in range(len(T)):
        if random.random()<pm:
            for j in range(n):
                num1,num2=integerToString()
            r=random.randint(0,len(strC)-2)
            if strC[r]=='0':
                strTemp=strC[:r]+'1'+strC[r+1:]
            else:
                strTemp=strC[:r]+'0'+strC[r+1:]
            x1=stringToInteger(strTemp[:len(strTemp)])
            x2 = stringToInteger(strTemp[len(strTemp):])
            R.append(Individual([x1,x2]))
        else:
            R.append(T[i])
    return R'''

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
def EMIA(Gmax,nd,na,nc):
    D=Initial(nd)
    for i in range(Gmax):
        print i
        DT=UpdateDomination(D,nd,na)
        A=ActiveSelection(DT,na)
        CT=RaClone(A,nc)
        Ct=CrossOver(CT,A,0.8)
        D=Ct+DT
    DT = UpdateDomination(D, nd, na)
    result=[]
    for i in range(len(DT)):
        if DT[i].bd==0:
            result.append(DT[i])
    x = []
    y = []
    for i in range(len(result)):
        x.append(result[i].f[0])
        y.append(result[i].f[1])
    pl.plot(x, y, '*')
    pl.xlabel('f1')
    pl.ylabel('f2')
    pl.show()


if __name__ == '__main__':
    EMIA(500,1200,400,1200)


