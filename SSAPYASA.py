from IPython.display import clear_output
import random
import time
from sage.all import *
p=2
m=1
q=p**m
F = GF(q,'a')
a=F.gen()
FL=F.list()
whatweight=0
dualflag=0
altref=1
ww2=1

def trail(t):
    r=0
    while t%2==0:
        r+=1
        t=t>>1
    return r



def weightenum(Cgen):
    if whatweight==0:
        if q==2:
            return weightenum1(Cgen)
        return weightenumq(Cgen)
    elif whatweight==1:
        return genweight(Cgen)
    elif whatweight==2:
        return smallword(Cgen)
    else:
        return biweight(Cgen)


    
def weightenum1(Cgen):
    n=Cgen.ncols()
    w=[0]*(n+1)
    w[0]+=1
    Cgen=Cgen.rref()
    while Cgen.nrows()!=0 and Cgen[-1]==vector([0]*n): 
        Cgen=Cgen[:-1]

    k=Cgen.nrows()

    c=[0]*n
    for s in range(1,2**(k)):
        i=trail(s)
        for l in range(n):
            c[l]=int(c[l])^int(Cgen[i][l]) #python
        k2=0
        for u in range(n):
            k2+=c[u]   
        w[k2]+=1
    return tuple(w)

def weightenumq(Cgen):
    n=Cgen.ncols()
    w=[0]*(n+1)
    w[0]+=1
    Cgen=Cgen.rref()
    while Cgen.nrows()!=0 and Cgen[-1]==vector([0]*n): 
        Cgen=Cgen[:-1]
    
    k=Cgen.nrows()
  
    for s in range(1,q**(k)):
        traq=[]
      
        while s!=0:
            traq+=[s%q]
            s=s//q
        
        traq+=(k-len(traq))*[0]
        c=[0]*n
       
        for u in range(k):
            for l in range(n):
                
                c[l]=(c[l]+FL[traq[u]]*Cgen[u][l])
        k2=0
        for u in range(n):
            if c[u]:
                k2+=1 
        w[k2]+=1
    return tuple(w)

def genweight(Cg):
    Cg=Cg.rref()
    w=[]
    while Cg.nrows()!=0 and Cg[-1]==vector([0]*n):
        Cg=Cg[:-1]
    
    n6=Cg.ncols()
    k6=Cg.nrows()
  
    c=zero_vector(n)
    cc=[]
    for s in range(1,2**(k6)):
        i=trail(s)
        for l in range(n):
            c[l]=int(c[l])^int(Cg[i][l]) #python
        cc+=[copy(c)]


    for k in range(1,k6+1):
        zv=transpose(matrix([0]*k))
        S=Subsets(range(len(cc)),k) 
        cand10=10000 
        for s in S:
            w6=0
            D=zero_matrix(0, n)
            for es in s:
                D=D.stack(cc[es])
            for j in range(n6):
                if D[:,j]!=zv:
                    w6+=1
            if w6<cand10:
                cand10=w6
        w+=[cand10]
    return tuple(w)
        
def smallword(Cg):
    Cg=Cg.rref()
    w=[]
    while Cg.nrows()!=0 and Cg[-1]==vector([0]*n): 
        Cg=Cg[:-1]
   
    n6=Cg.ncols()
    k6=Cg.nrows()
    if k6==0:
        return 0
    c=zero_vector(n)
    cc=[]
    for s in range(1,2**(k6)):
        i=trail(s)
        for l in range(n):
            c[l]=int(c[l])^int(Cg[i][l]) #python
        k2=0
        for u in range(n):
            k2+=c[u] 
        if k2 not in cc:
            cc+=[k2]
    cc.sort()
    cc=cc[:ww2]
    return cc[-1]

def biweight(Cg):
    Cg=Cg.rref()
    n=Cg.ncols()
    w=[0]*(n+1)
    w[0]+=1
    while Cg.nrows()!=0 and Cg[-1]==vector([0]*n):
        Cg=Cg[:-1]
  
    n6=Cg.ncols()
    k6=Cg.nrows()
   
    c=zero_vector(n)
    cc=[]
    for s in range(1,2**(k6)):
        i=trail(s)
        for l in range(n):
            c[l]=int(c[l])^int(Cg[i][l]) #python
        cc+=[copy(c)]
        k2=0
        for u in range(n):
            k2+=c[u]  
        w[k2]+=1
        


    for k in range(1,k6+1):
        zv=transpose(matrix([0]*k))
        S=Subsets(range(len(cc)),k)
        cand10=10000 
        for s in S:
            w6=0
            D=zero_matrix(0, n)
            for es in s:
                D=D.stack(cc[es])
            for j in range(n6):
                if D[:,j]!=zv:
                    w6+=1
            if w6<cand10:
                cand10=w6
        w+=[cand10]
    return tuple(w)
        

def dsf(B,A=None):
    return dsf2(B,A)


def dsf2(B,A=None):
    #print("B:")
    #print(B)
    #print("Bend")
    '''
    B=B.rref()
    global n
    n=B.ncols()
    while B[-1]==vector([0]*n):
        B=B[:-1]
    '''
    B2=B.stack(zero_matrix(n-B.nrows(), n))

    B2=B2.augment(identity_matrix(F,n))
    M2=B2.rref() 
   
    R=M2[:,:n]
    T=M2[:,n:]
    while R[-1]==vector([0]*n):
        R=R[:-1]
    kd2=R.nrows()
    i=0
    while i<kd2 and kd2<n:
        if R[i][i]==0:
            R=R[:i].stack(vector([0]*n)).stack(R[i:])
            T=T[:i].stack(T[kd2]).stack(T[i:])
            T=T[:kd2+1].stack(T[kd2+2:])
            kd2+=1
        i+=1
    R=R.stack(zero_matrix(n-kd2, n))

    
    if A is not None:
        L=T*A
    
    else:
        L=None
    
    return L,R,T


def dsf3(B,A=None):
    
    if A==None:
        A=identity_matrix(F,n,n)
    L=copy(A)
    R=B.stack(zero_matrix(F,n-B.nrows(), n))
    R=R.rref()
    R=copy(R)
    for k in range(n):

        if R[k,k]==0:
            for i in range(n):
                if i==k:
                    continue
                if R[k,k]==0 and (i>=k or R[i,i]==0) and R[i,k]!=0:
                    L[k]=L[k]+L[i]
                    R[k]=R[k]+R[i]
                    break
            if R[k,k]==0:
                continue
        L[k]=(1/R[k,k])*L[k]
        R[k]=(1/R[k,k])*R[k]
        for i in range(n):
            if i==k:
                continue
            if R[i,k]==0:
                continue
            L[i]=L[i]-R[i,k]*L[k]
            R[i]=R[i]-R[i,k]*R[k]
    return L,R



def hully(M):
    
    A=zero_matrix(F,n,n) 
    B=copy(M) 
    wherepc=[]
 
    for k in range(n):
        if M[k,k]==0:
            wherepc+=[k]
            A[k,k]=1 
            B[k,k]=1
            for j in range(k):
                A[k,j]=-M[j,k]
                B[k,j]=-M[j,k]
  
    L,R,T=dsf(B,A)
    l=-1
  
    H=zero_matrix(F,0,n) 
    for k in range(n):
        if R[k,k]==0:
            H=H.stack(vector([0]*n))
            l+=1
            H[l,k]=1
            for j in range(k):
                H[l,j]=-R[j,k]

    return H,L,B,R,wherepc,A,T


def phull(H,L,r):
    Hp=copy(H)
    l=H.nrows()
    for k in reversed(range(l)):
        if H[k,r]:
            for i in range(k):
                if H[i,r]:
                    
                    Hp[i]=H[i]+H[k]
                   
            Hp=Hp.delete_rows([k])
            return Hp
   
    Hp=Hp.stack(L[r])
    

    return Hp



def evfilter(FF,C,flaggy=0):
 
    H,L=hully(dsf(C)[1])
    up=weightenum(H)
    if FF[0]!=up:
        if flaggy==1:
            print("Hull wrong")
        return False
    l=len(FF[1])
    vp=[0]*n
    for k in range(n):
        vp[k]=weightenum(phull(H,L,k))
        for i in range(l):
         
            if FF[1][i][0]==False and vp[k]==FF[1][i][1]:
                print("forbidden punc")
                return False
    for i in range(l):
        if FF[1][i][0]:
            r=False
            for k in range(n):
               
                if vp[k]==FF[1][i][1]:
                    r=True
                    break
            if not r:
                if flaggy==1:
                    print("missing punc")
                return False
    return True

def buildfilter(C,K=100000):
    H,L=hully(dsf(C)[1])
    u=weightenum(H) 
    Vy=[]
    for k in range(n):
        Vy.append(weightenum(phull(H,L,k)))
    Cl=[]
    V=Set(Vy)
    for k in range(K):
        Cidk = random_matrix(F,kd,n)
        Cl.append(Cidk)
        Vtemp=[]
        for j in range(n):
            Vtemp.append(weightenum(phull(H,L,j)))
        V=V.union(Set(Vtemp)) 
    FF=[u,[]]
    while len(V)>0:
        for C in Cl:
            if not evfilter(FF,C):
                Cl.remove(C)
               
        if Cl==[]:
            print("None left")
            break
        epsar=[]
        for v in V:
           
            b=False
            if v in Vy:
                b=True
            eps=0
            for C in Cl:
                if evfilter([u,[[b,v]]],C):
                    eps+=1
            epsar.append([eps,v,b])
        epsar=sorted(epsar, key=lambda student: student[0]) 
        V=set(V)
        V.remove(epsar[0][1])
        V=Set(V)
        FF[1].append([epsar[0][2],epsar[0][1]])
    return FF


def SSAprep(C):
   
    H,L,M2,D2,wherepc2,A2,T=hully(C)
    
    u=weightenum(H) 
    Vy=[]
    it=1

    Cd=dsf(dualc(C))[1] 
    _,Ld,_,_,_,_,Td=hully(Cd)
    Cdp=[]
    
    for k in range(n):
        flag10=0
       
        
        if altref==0:
            w=weightenum(phull(H,L,k))+weightenum(phull(H,Ld,k))#dull
        else:
           
            wi=phull(H,L,k) 
           
            w=weightenum(wi)
           
        if [w,0,0] in Vy:
            
            Vy[Vy.index([w,0,0])]=[w,it,0]
            Vy.append([w,it,0])
            it+=1
        else:
            for j in range(1,it):
                if [w,j,0] in Vy: 
                    Vy.append([w,j,0])
                    flag10=1
                    break
            if flag10==0:
                Vy.append([w,0,0])

        Cpc=puncc(C,k)
        Cdpc=puncc(Cd,k)
        Cpcd=dualc(dsf(Cpc)[1])
        Cdpcd=dualc(dsf(Cdpc)[1])
        Cdp+=[[Cpcd,Cdpcd]]
  
    Vy.append([u,-1]) 
   
    return Vy,H,T,Td,Cd,Cdp,it,M2,D2,wherepc2,A2
    
def SSA(C,prepss,debugC):
    
    if rank(C)!=rank(debugC):
        return "dim mismatch" 
    preps=deepcopy(prepss)
    C=dsf(C)[1]
    debugC=dsf(debugC)[1]
    
    prep,sH,sL,sLd,sCd,sCdp,_,sM2,sD2,swherepc2,sA2=preps #TODO remove unnecessary
    Vy,H,L,Ld,Cd,Cdp,it,M2,D2,wherepc2,A2=SSAprep(C)

    
    
    if Vy[-1]!=prep[-1]:
        return "hull wrong"
    
    prep=prep[:-1]
    Vy=Vy[:-1]
   
    timer=1
   
    perm=[-1]*(n)
    ird=[]

    brtlst=[]
  
    prepslice=[x[0] for x in prep]
    Vslice=[x[0] for x in Vy]
    for k in range(n):
        if Vy[k][0] not in prepslice: 
          
            print("failed at",k)
            return "no presence, filtered"
        if prepslice[k] not in Vslice:
          
            return "no presence 222"
        if Vy[k][1]==0:
            if Vy[k] not in prep:
                print("failed at",k)
                return "jhkasasfh"
            perm[k]=prep.index(Vy[k])
            ird+=[k]
            
    bruteflag=[]
    yaf2=[]
    yaf3=[] 
    Vimprint=[] 
    permimprint=[]
    brtlstimprint=[]
    prepimprint=[]
    irdimprint=[]
    Vtemp="def" 
    Vtemp2=[]
    sVtemp="ghi" 
    sVtemp2=[] 
    
    YAF=0
    yaf22=0
    yaf55=0
    
    yafref=0 
    yafref2=0
    ref1=0
    ref2=0
    
    flag90=0
    dif=1
    
    while True:
        #clear_output()
        #print("irdyy:",ird,"\n",[x[1] for x in Vy])
        
        if yaf22==1:
            #print("brutyyy,max:",bruteflag,yaf3,"\n",perm,"\n",debugperm2)
            yaf22=0
            timer=1
            for i in reversed(range(1,len(bruteflag))): 
                if bruteflag[i]>=yaf3[i]: 
                    bruteflag=bruteflag[:-1]

                    bruteflag[-1]+=1
                    
                    yaf3=yaf3[:-1]
                    YAF=1
                   
                    Vimprint=Vimprint[:-1]
                   
                    prepimprint=prepimprint[:-1]
                    
                    irdimprint=irdimprint[:-1]
                    
                    
                    permimprint=permimprint[:-1] 

                    brtlstimprint=brtlstimprint[:-1]
                   

                    Vtemp2=Vtemp2[:-1]
                    sVtemp2=sVtemp2[:-1]

        
            if YAF==1:
                print("aww")
                yafref2=1
                if bruteflag!=[]:
                    if bruteflag[0]>=yaf3[0]:
                        print("breaks for large q sometimes")
                        return "out of guesses"
                else:
                    return "exhausted"
                YAF=0
                perm=copy(permimprint[-1])            
                prep=copy(prepimprint[-1]) 
                Vy=copy(Vimprint[-1])
                ird=copy(irdimprint[-1])
                brtlst=copy(brtlstimprint[-1]) 
        
                perm[Vtemp2[-1][0]]=sVtemp2[-1][bruteflag[-1]] 
                prep[sVtemp2[-1][bruteflag[-1]]][1]=0
                Vy[Vtemp2[-1][0]][1]=0
           
               
            else: 
                print("decided: else, going deeper")
                slice1=[y for (x,y,z) in Vy]

                county=[slice1.count(itt) for itt in range(1,it)] 
                for iter49 in range(len(county)):
                    if county[iter49]==0:
                        county[iter49]=10000 
        
                it6=county.index(min(county))+1 
                
                
                tempadd=[i for i, x in enumerate(Vy) if x[1] == it6]
                if tempadd==[]:
                    return "no more"
               
                
                cand3=findmatch(Vy,tempadd[0],prep)
               
                if cand3=="nah":
                    return "no more"
                if cand3=="nah2":
                     
                    return "no more2"
                sVtemp2+=[cand3]
                Vtemp2+=[tempadd]
            
                if len(Vtemp2[-1])!=len(sVtemp2[-1]):

                    return "temps not same len1"
                
                bruteflag+=[0] 
                
                yaf3+=[len(Vtemp2[-1])]
    
                permimprint+=[copy(perm)]
                
                brtlstimprint+=[brtlst] 
                prepimprint+=[prep]
                irdimprint+=[ird]
                Vimprint+=[Vy]
               
                
                perm[Vtemp2[-1][0]]=sVtemp2[-1][0] 
                Vy[Vtemp2[-1][0]][1]=0
                prep[sVtemp2[-1][0]][1]=0
                brtlst+=[Vtemp2[-1][0]]
              
        
        irdtemp=[x for x in ird if randrange(4)==0]+brtlst
        irdtemp2=[perm[u] for u in irdtemp]
        varra=list(range(1,it+1))

        for var24 in varra:
            if randrange(4):
                continue
            adder=[i for i,x in enumerate(Vy) if x[1] == var24]
            if adder==[]:
                continue
            irdtemp+=adder

            cand=findmatch(Vy,adder[0],prep)
            if cand=="nah":
                continue
            if cand=="nah2":
                yaf22=1
                YAF=1
                bruteflag[-1]+=1
                yaf55=1
                break
            irdtemp2+=cand
      
        if yaf55==1:
            yaf55=0
            timer+=1
            continue
        if irdtemp==[]:
           
            timer+=1
            continue
        
        for it3 in range(1,it):
            if YAF==1:

                break
           
            Vtemp=[i for i, x in enumerate(Vy) if x[1] == it3]
            if Vtemp==[]:
                continue
            
            sVtemp=findmatch(Vy,Vtemp[0],prep)
            
            if sVtemp=="nah":
               
                continue
            if sVtemp=="nah2":
               
                
                yaf22=1
                YAF=1
                bruteflag[-1]+=1
                break
                
            if Vtemp[0] in irdtemp:
                for u in Vtemp:
                    irdtemp.remove(u)
                for v in sVtemp:
                    irdtemp2.remove(v)
                flag90=1
          
            if len(Vtemp)!=len(sVtemp):
                
                if bruteflag!=[]:
                    YAF=1
                    yaf22=1
                    bruteflag[-1]+=1
                    continue
                
                return "temps not same len2"
          
         
            V1,it1=refine(irdtemp,Vtemp,C,Cd,Cdp,it,L,Ld,M2,D2,wherepc2,A2)
            V2,sit1=refine(irdtemp2,sVtemp,debugC,sCd,sCdp,it,sL,sLd,sM2,sD2,swherepc2,sA2)
            if it1!=sit1:
                if bruteflag==[]:
                    return "refine fail"
              
                yaf22=1
               
                YAF=1
                bruteflag[-1]+=1
                break
                
            if it1==it or it1>it+1: 
                yafref=1
            
            it=it1+1
        
            if flag90==1:
                flag90=0
                irdtemp+=Vtemp
                irdtemp2+=sVtemp
                
            for k in range(len(Vtemp)):
        
                V2slice=[[x[0],x[2]] for x in V2]
                V1slice=[[x[0],x[2]] for x in V1]
                if V1slice[k] not in V2slice:
                    if bruteflag==[]:
    
                        return "not equiv after refine"
        
                    yaf22=1
                  
                    YAF=1
                    bruteflag[-1]+=1
                    break
              
            
                if V1[k][1]==0:
                    yafref=1
                        
                    perm[Vtemp[k]]=sVtemp[V2.index(V1[k])]
                   
                    ird+=[Vtemp[k]]
                  
                Vy[Vtemp[k]][1]=V1[k][1]
                prep[sVtemp[k]][1]=V2[k][1]
                Vy[Vtemp[k]][0]=V1[k][0]
                prep[sVtemp[k]][0]=V2[k][0]
                Vy[Vtemp[k]][2]=dif
                prep[sVtemp[k]][2]=dif
            dif+=1
                
        
        timer+=1
        if yafref==1:
            yafref=0
            ref1+=1
      
        ref2+=1
        
        if -1 not in perm:
            pmre=[i+1 for i in perm]
            pemr = Permutation(pmre).inverse().to_matrix()
            if dsf(C*pemr)[1]==debugC:
              
                if yafref2==1:
                    return perm,0,0
                return perm,ref1,ref2
            
            YAF=1
            bruteflag[-1]+=1
            yaf22=1
           


        if bruteflag==[]:
            permeq=copy(perm)
            
            for it9 in range(1,it):
                Vtemp9=[i for i, x in enumerate(Vy) if x[1] == it9]
                if Vtemp9==[]:
                    continue
                sVtemp9=findmatch(Vy,Vtemp9[0],prep)
                if sVtemp9=="nah2":
                  
                    yaf22=1
                    YAF=1
                    bruteflag[-1]+=1
                    yaf55=1
                    break
                
                for h in range(len(Vtemp9)):
                    permeq[Vtemp9[h]]=sVtemp9[h]
            if yaf55==1:
                yaf55=0
                timer+=1
                continue
          
            pmre=[i+1 for i in permeq]
            pemr = Permutation(pmre).inverse().to_matrix()
            if dsf(C*pemr)[1]==debugC:
                print("refseq:",ref1,ref2)
                if yafref2==1:
                    return permeq,0,0
                return permeq,ref1,ref2
        

        
        if timer>=50:
            yaf22=1
            print("time reached; guessing")

  
def refine(pnc,Vtemp,DC,DDC,Cm,it7,L,Ld,M2,D2,wherepc2,A):
    if altref==0:
        return oldrefine(pnc,Vtemp,DC,DDC,Cm,it7)
    return newrefine(pnc,Vtemp,DC,DDC,Cm,it7,L,Ld,M2,D2,wherepc2,A)

def oldrefine(pnc,Vtemp,DC,DDC,Cm,it7):
    Vrf=[]
    it=it7
    if dualflag==0:
        for t in pnc:
            DC=puncc(DC,t)
    else:
        for t in pnc:
            DC=puncc(DC,t)
            DDC=puncc(DDC,t)
    for i in Vtemp: 
        w11=[]
        inter1=puncc(DC,i)
        if dualflag==1:
            inter2=puncc(DDC,i) 

            Cm1=Cm[i][0]
            Cm2=Cm[i][1]
            for t in pnc:
                Cm1=puncc(Cm1,t)
                Cm2=puncc(Cm2,t)
        
            for mat in [inter1,inter2,Cm1,Cm2]:
                if mat==zero_matrix(F,kd,n) or mat==zero_matrix(F,n,n):
                    debby=matrix(F,[])
                else:
                    debby=hully(dsf(mat)[1])[0]
                w11+=weightenum(debby)
        elif dualflag==2:
            inter2=puncc(DDC,i) 
            for mat in [inter1,inter2]:
                if mat==zero_matrix(F,kd,n) or mat==zero_matrix(F,n,n):
                    debby=matrix(F,[])
                else:
                    debby=hully(dsf(mat)[1])[0]
                w11+=weightenum(debby)
        else:
            if inter1==zero_matrix(F,kd,n) or inter1==zero_matrix(F,n,n):
                debby=matrix(F,[])
            else:
                debby=hully(dsf(inter1)[1])[0]
            w11+=weightenum(debby)
        flag10=0
        
        if [w11,0,1] in Vrf: 
            Vrf[Vrf.index([w11,0,1])]=[w11,it,2]
            Vrf.append([w11,it,2])
            it+=1
        else:
            Vrf2=[[x[0],x[1]] for x in Vrf]
            for j in range(it7,it):
                if [w11,j] in Vrf2: 
                    wr=Vrf2.index([w11,j])
                    new=Vrf[wr][2]+1
                    for h in Vrf:
                        if h[1]==j:
                            h[2]=new
                    Vrf.append([w11,j,new])
                    flag10=1
                    break
            if flag10==0:
                Vrf.append([w11,0,1])
    
    return Vrf,it

def newrefine(pnc,Vtemp,DC,DDC,Cm,it7,L,Ld,M,D,wp,A):
    Vrf=[]
    it=it7

    MP=copy(M)

    MP2=copy(M)
  
    for j in range(n):
        if j in wp:
            for t in pnc:
                MP2[j,t]=0
            continue
        for t in pnc:
            MP[j,t]=0
  
    X=MP-M

    X2=MP2-M
    SP,DP,_=dsf(D+L*X,L)

    LP=SP*A
    
    SP2,DP2,_=dsf(D+L*X,L)
    LP2=SP2*DC
   

    l=-1
    l2=-1
    BP=zero_matrix(F,0,n) 
    BP2=zero_matrix(F,0,n)
    for k in range(n):
        if DP[k,k]==0:
            BP=BP.stack(vector([0]*n))
            l+=1
            BP[l,k]=1
            for j in range(k):
                BP[l,j]=-DP[j,k]
        if DP2[k,k]==0:
            BP2=BP2.stack(vector([0]*n))
            l2+=1
            BP2[l,k]=1
            for j in range(k):
                BP2[l,j]=-DP2[j,k]
   
    for t in pnc:
        BP=puncc(BP,t)
        BP2=puncc(BP2,t)
       
    for i in Vtemp:
        
        predeb=phull(BP,LP,i)
        #predeb2=phull(BP2,LP2,i)
      
        w11=weightenum(predeb)#+weightenum(predeb2)
        flag10=0

        
        if [w11,0,1] in Vrf:
            Vrf[Vrf.index([w11,0,1])]=[w11,it,2]
            Vrf.append([w11,it,2])
            it+=1
        else:
            Vrf2=[[x[0],x[1]] for x in Vrf]
            for j in range(it7,it):
                if [w11,j] in Vrf2: 
                    wr=Vrf2.index([w11,j])
                    new=Vrf[wr][2]+1
                    for h in Vrf:
                        if h[1]==j:
                            h[2]=new
                    Vrf.append([w11,j,new])
                    flag10=1
                    break
            if flag10==0:
                Vrf.append([w11,0,1])
    
    return Vrf,it


    
def findmatch(V1,some,V2):
  
    ind2=-1
   
    currph=V1[some][0]
    currz=V1[some][2]
    
    for (x,y,z) in V2:
        if x==currph and y!=0 and z==currz: 
            ind2=y
            break
    if ind2==-1:
      
        return "nah2"
   
    return [i for i,(x,y,z) in enumerate(V2) if y==ind2]    


def squarec(C):
    k=C.nrows()
    print("k:",k)
    PowC=zero_matrix(F,0,n)
    for i in range(k):
        for j in range(i+1):
           
            temprow=zero_vector(F,n)
            for l in range(n):
                temprow[l]=C[i,l]*C[j,l]
             
            PowC=PowC.stack(temprow)
           
    PowC=PowC.rref()
    while PowC.nrows()!=0 and PowC[-1]==vector([0]*n):
        PowC=PowC[:-1]
   
    return PowC

def puncc(C,pnc):
    C6=copy(C)
    
    k8=C6.nrows()

    C6[:,pnc]=transpose(matrix([0]*k8))
   
    return C6



def extc(C):
   
    return "nah"

def dualc(R):
    
    l=-1
    H=zero_matrix(F,0,n)
    for k in range(n):
        if R[k,k]==0:
            
            H=H.stack(vector([0]*n))
            l+=1
            H[l,k]=1
            for j in range(k):
                H[l,j]=-R[j,k]
    return H



    
def iterator(C1,n,k3,blac=[]):
    prip=SSAprep(dsf(C1)[1])
    subs1=list(range(n))
    S=Subsets(subs1,k3)
    count3=0
    count4=0
    bla=[]
    print("beforeqbin:n,k:",n,k3)
    print("qbin:",gaussian_binomial(n,k3,2))
    for s in S:
        #print(s)
        A=zero_matrix(F,n,n)
        for k in s:
            A[k,k]=1
        #print("new")
        g2=[]
        l2=[]
        l3=0
        for k in s:
            g=list(range(k+1,n))
            g2+=[list(set(g) - set(s))]
            l2+=[len(g2[-1])]
            l3+=l2[-1]
            #print("g,l:\n",g2,l2,l3)
        print("maxx:",2**(l3))
        for h in range(2**(l3)):
            #clear_output()
            #print("h",h)
            h+=0 #how to convert int to sageint
            bi1=h.digits(2)
            #print("bi1:",bi1)
            bi1.reverse()
            bi=[0]*(l3-len(bi1))+bi1# or just keep it reversed b+0
            #print("bi:",bi)
            ctr1=0
            for j in range(len(s)):
                for d in g2[j]:
                    A[s[j],d]=bi[ctr1]
                    ctr1+=1
            if A in blac:
                count3+=1
                continue
            A2=copy(A) #deepcopy
            prap=copy(prip) #deepcopy
            #print("A2:",A2)
            ess=SSA(A2,prap,C1)
            #print("trying:",A2)
            if type(ess)==list:
                #print("check:\n",ess,ess==True)
                #print("success:",C1,"\n\n",A2)
                count4+=1
                bla+=[A2]
            else:
                count3+=1
        
    print("count1,2:",count4,count3)
    return count4,bla
    
def iterator2(n,k2): #dont use
    exit()
    subs1=list(range(n))
    S=Subsets(subs1,k2)
    bla=[]
    amount=[]
    print("beforefakeqbin:n,k:",n,k2)
    print("fakeqbin:",gaussian_binomial(n,k2,2))
    for s in S:
        #print(s)
        A=zero_matrix(F,n,n)
        for k in s:
            A[k,k]=1
        #print("new")
        g2=[]
        l2=[]
        l3=0
        for k in s:
            g=list(range(k+1,n))
            g2+=[list(set(g) - set(s))]
            l2+=[len(g2[-1])]
            l3+=l2[-1]
            #print("g,l:\n",g2,l2,l3)
        print("fakemaxx:",2**(l3))
        for h in range(2**(l3)):
            #clear_output()
            #print("h",h)
            h+=0 #how to convert int to sageint
            bi1=h.digits(2)
            #print("bi1:",bi1)
            bi1.reverse()
            bi=[0]*(l3-len(bi1))+bi1# or just keep it reversed b+0
            #print("bi:",bi)
            ctr1=0
            for j in range(len(s)):
                for d in g2[j]:
                    A[s[j],d]=bi[ctr1]
                    ctr1+=1
            if A in bla:
                continue
            A2=copy(A) #deepcopy
            am1,bla1=iterator(A2,n,k2,bla)
            amount+=[am1]
            bla+=bla1
            #prap=deepcopy(prip)
            #print("A2:",A2)
            #ess=SSA(A2,prap,C1)
            #print("trying:",A2)
            #if type(ess)==list:
                #print("check:\n",ess,ess==True)
                #print("success:",C1,"\n\n",A2)
             #   count4+=1
             #   bla+=[A2]
            #else:
             #   count3+=1
        
    print(amount)
    return amount

def iterator3(n,k4): #GF? with duplicates as to not falsify distribution
    subs1=list(range(n))
    S=Subsets(subs1,k4)
    count3=0
    count4=0
    bla=[]
    amount=[]
    yaf24=0
    print("qbin:",gaussian_binomial(n,k4,2))
    for s in S:
        #print(s)
        A=zero_matrix(F,n,n)
        for k in s:
            A[k,k]=1
        #print("new")
        g2=[]
        l2=[]
        l3=0
        for k in s:
            g=list(range(k+1,n))
            g2+=[list(set(g) - set(s))]
            l2+=[len(g2[-1])]
            l3+=l2[-1]
            #print("g,l:\n",g2,l2,l3)
        #print("maxx:",2**(l3))
        for h in range(2**(l3)):
            #if h%1000==0:
            #    print(h)
            #clear_output()
            #print("h",h)
            h+=0 #how to convert int to sageint
            bi1=h.digits(2)
            #print("bi1:",bi1)
            bi1.reverse()
            bi=[0]*(l3-len(bi1))+bi1# or just keep it reversed b+0
            #print("bi:",bi)
            ctr1=0
            for j in range(len(s)): #or len(s)=k4
                for d in g2[j]:
                    A[s[j],d]=bi[ctr1]
                    ctr1+=1
            prip=SSAprep(dsf(A)[1])
            for h in bla:
                A2=copy(A) #deepcopy
                #prap=deepcopy(prip)
                ess=SSA(h,prip,A2)
                if type(ess)==list:
                    amount[bla.index(h)]+=1
                    yaf24=1
                    break
            if yaf24==1:
                yaf24=0
                continue
            #print("newclass")
            bla+=[copy(A)]#deepcopy
            amount+=[1]
            #if len(amount)==7:
            #    print("A:",bla[6])
            #    print(amount)
            #    exit()
            
            #A2=deepcopy(A)
            #prap=deepcopy(prip)
            #print("A2:",A2)
            #ess=SSA(A2,prap,C1)
            #print("trying:",A2)
            
        
    print("n,k,classes",n,k4,amount)
    return bla,amount

def goppa3(n,k4,bla): #defunct 
    exit()
    subs1=list(range(n))
    S=Subsets(subs1,k4)
    count3=0
    count4=0
    bla=[]
    amount=[0]*len(bla)
    yaf24=0
    print("qbin:",gaussian_binomial(n,k4,2))
    for s in S:
        #print(s)
        A=zero_matrix(F,n,n)
        for k in s:
            A[k,k]=1
        #print("new")
        g2=[]
        l2=[]
        l3=0
        for k in s:
            g=list(range(k+1,n))
            g2+=[list(set(g) - set(s))]
            l2+=[len(g2[-1])]
            l3+=l2[-1]
            #print("g,l:\n",g2,l2,l3)
        #print("maxx:",2**(l3))
        for h in range(2**(l3)):
            #if h%1000==0:
            #    print(h)
            #clear_output()
            #print("h",h)
            h+=0 #how to convert int to sageint
            bi1=h.digits(2)
            #print("bi1:",bi1)
            bi1.reverse()
            bi=[0]*(l3-len(bi1))+bi1# or just keep it reversed b+0
            #print("bi:",bi)
            ctr1=0
            for j in range(len(s)): #or len(s)=k4
                for d in g2[j]:
                    A[s[j],d]=bi[ctr1]
                    ctr1+=1
            prip=SSAprep(dsf(A)[1])
            for h in bla:
                A2=copy(A)#deepcopy
                #prap=deepcopy(prip)
                ess=SSA(h,prip,A2)
                if type(ess)==list:
                    amount[bla.index(h)]+=1
                    yaf24=1
                    break
            #if yaf24==1:
            #    yaf24=0
            #    continue
            #print("newclass")
            #bla+=[deepcopy(A)]
            #amount+=[1]
            #if len(amount)==7:
            #    print("A:",bla[6])
            #    print(amount)
            #    exit()
            
            #A2=deepcopy(A)
            #prap=deepcopy(prip)
            #print("A2:",A2)
            #ess=SSA(A2,prap,C1)
            #print("trying:",A2)
            
        
    print("goppa:n,k,classes",n,k4,amount)
    return bla,amount


    
def splitter(C):
    print("splitlimit 10 cause runtime")
    H,L=hully(dsf(C)[1])
    Vy=[]
    it=1
    #print("H:",H)
    for k in range(n):
        flag10=0
        w=weightenum(phull(H,L,k))
        #print("w:",w,phull(H,L,k))
        if [w,0] in Vy:
            
            Vy[Vy.index([w,0])]=[w,it]
            Vy.append([w,it])
            it+=1
        else:
            for j in range(1,it):
                if [w,j] in Vy: #cleanup better elif
                    Vy.append([w,j])
                    flag10=1
                    break
            if flag10==0:
                Vy.append([w,0])
    ird7=[Vy.index(x) for x in Vy if x[1]==0]
    print("before:",[x[1] for x in Vy])
    print("ird7 1:",ird7)
    VyC=deepcopy(Vy)
    timer=0
    flag90=0
    while(True):
        #print("ird at time",ird7,timer)
        #print([x[1] for x in Vy])
        irdtemp=[x for x in ird7 if randrange(4)==0]
        varra=list(range(1,it+1))
        for var24 in varra:
            if randrange(4):
                continue
            adder=[i for i, x in enumerate(Vy) if x[1] == var24]
            if adder==[]:
                continue
            irdtemp+=adder
        if irdtemp==[]:
            continue

        for it3 in range(1,it):
            
            '''
            irdt=deepcopy(ird7)
            #shuffle(irdt)
            #irdtemp=irdt[:(len(irdt)//2)]#TODO all with 1/2
            irdtemp=[x for x in irdt if randrange(2)==1]
            varra=list(range(1,it+1))
            #print("varra:",varra)
            varra.remove(it3)
             3x
            while len(irdtemp)<timer:
                if varra==[]:
                    #print("whatever go my scarab")
                    #yaf22=1
                    break
                #ggprint("whiling")
                var24=random.choice(varra)
                #ggprint("g is",var24)
                varra.remove(var24)
                adder=[i for i, x in enumerate(Vy) if x[1] == var24]
                #ggprint("addy:",adder)
                if adder==[]:
                    continue
                irdtemp+=adder
             3x
            for var24 in varra:
                if randrange(2):
                    continue
                adder=[i for i, x in enumerate(Vy) if x[1] == var24]
                if adder==[]:
                    continue
                irdtemp+=adder
            #print("ord::::",irdtemp)
            if irdtemp==[]:
                #print("AAAAAAAAAAAA",len(irdtemp)<timer)
                #itprint("ird is:",ird)
                #print("irdnahhh")
                #print(C)
                break
            '''

                
            Vtemp=[i for i, x in enumerate(Vy) if x[1] == it3]
            if Vtemp==[]:
                continue
            if Vtemp[0] in irdtemp:
                for u in Vtemp:
                    irdtemp.remove(u)
                flag90=1
            V1,it1=refine(irdtemp,H,L,Vtemp,C,it)
            it=it1+1
            for k in range(len(Vtemp)):
                if V1[k][1]==0:
                     ird7+=[Vtemp[k]]
                Vy[Vtemp[k]][1]=V1[k][1]
                Vy[Vtemp[k]][0]=V1[k][0]
            if flag90==1:
                flag90=0
                irdtemp+=Vtemp
        #print("Vyafter:",Vy)
        timer+=1
        if len(ird7)==n:
            print("early finish with timer:",timer)
            break
        if timer>=5:
            break
    print("ird7 2 :",ird7)
    #print("broodfoaß")
    print("after:",[x[1] for x in Vy])
    
    '''
    varra=list(range(1,it+1))
    for i in range(1,it+1):
        print("brute card",i)
        S2=Subsets(varra,i)
        for it3 in range(1,it):
            for s in S2:
                if it3 in s:
                    continue
                Vtemp=[i for i, x in enumerate(Vy) if x[1] == it3]
                if Vtemp==[]:
                    break
                irdtemp7=[]
                for z in s:
                    irdtemp7+=[i for i, x in enumerate(VyC) if x[1] == z]


                    V7,it1=refine(irdtemp7,H,L,Vtemp,C,it)
                    it=it1+1
                    for k in range(len(Vtemp)):
                        if V7[k][1]==0:
                            Vy[Vtemp[k]][1]=0
                            ird7+=[Vtemp[k]]
    '''


            
    #print("finished:",Vy,"\n",ird7)
        
    '''
    H,L=hully(dsf(C)[1])
    V=list(range(n))
    pncz=[]
    V6=refine(pncz,H,L,V,C)'''
    #print(V6)

#def mGRS(q,n,k,alph,bet): #for now q=2 TODO F or GF(q) or GF(2)
#    A=zero_matrix(GF(q),k,n)
#    for i in range(n):
#        for j in range(k):
#            A[j,i]=bet[i]*pow(alph[i],j,q)
#    return A


#FF.<X> = FF[]
#print(FF)
#print(PFF)
#print(FF.gen()^8)
            
def mGoppa(g,L):
    #print("g:",g)
    D=len(g)
    C=zero_matrix(F,mg*(D-1), 0)
    for i in range(len(L)):
        G=0
        for j in range(D):
            G+=(L[i]**j)*g[j]
        #print("G:",G)
        if G==0:
            return "no can do"
        aux=[]
        for h in range(D-1):
            #print("beforeaux:",(L[i]^h)/G)
            aux+=((L[i]**h)/G).list()
        #print(aux)
        C=C.augment(vector(aux))
        #print("updatedC:",C)
    return C

def expand(C):
    n=C.ncols()
    k=C.nrows()
    C2=matrix(GF(p),m*k,0)
    print("mk:",m*k)
    for i in range(n):
        aux=[]
        for j in range(k):
            aux+=C[j,i].list()
        C2=C2.augment(vector(aux))
    return C2

def closure(C):
    n=C.ncols()
    k=C.nrows()
    C2=matrix(GF(p),k,0)
    
    for i in range(n):
        for j in FL:
            aux=[]
            for t in range(k):
                aux+=[j*C[t,i]]
            C2=C2.augment(vector(aux))
    return C2


def mGRS(L,bet):
    #global n,kd
    #kd=k6
    #n=len(L)
    C=zero_matrix(F,kd, 0)
    for i in range(n):
        aux=[]
        for h in range(kd):
            #print("ggg:",bet[i]*(L[i]^h))
            aux+=[bet[i]*(L[i]**h)]
            #print("newaux:",aux)
        #print("aux:",aux)
        C=C.augment(vector(aux))
        #print("Cnow:",C)
    return C


    #lup
#mGoppa([1,a,1],[0,1,a,a+1])
#print("FFFF:",FF.list())

#PFF.<x> = PolynomialRing(FF)
#print("grss:",mGRS([FF(1),a,a^2,a^3,a^4,a^5,a^6],[1,1,1,1,1,1,1],2))
#print(FF.list())
def GRSiterator():
    #print("dontusecurrently")
    #global kd
    #F1=F.list() #FF
    S1=Subsets(range(q),n)
    bettemp=[1]*n
    amount=[]
    bla=[]
    yaf24=0
    h0=0
    h1=0
    print(binomial(q,n))
    hk=[0]*kd
    #for s in S1:
    for _ in range(10000):
        s=S1.random_element()
        bettemp=[]
        for i in range(n):
            bettemp+=[randrange(1,n)]
        #print(len(bettemp))
        
        L1=[]
        for z in s:
            L1+=[FL[z]]
        A=mGRS(L1,bettemp)
        hd=hully(dsf(A)[1])[0].nrows()
        hk[hd]+=1
        if _==50:
            print("check")
    print(hk)
'''
        A=mGRS(L1,bettemp,kd).rref()
        while A[-1]==vector([0]*n): #dont know where to place
            A=A[:-1]
        A=copy(A)
        kd=rank(A)
        prip=SSAprep(dsf(A)[1])
        for h in bla:
            A2=copy(A)#deepcopy
            #prap=deepcopy(prip)
            ess=SSA(h,prip,A2)
            if type(ess)==list:
                amount[bla.index(h)]+=1
                yaf24=1
                break
        if yaf24==1:
            yaf24=0
            continue
        bla+=[copy(A)]#deepcopy
        amount+=[1]
    print("blaamount:",bla,amount)
    #print("nk",n,kd)
    amount2=[]
    for A in bla:
        A2=copy(A)#deepcopy
        amount2+=[iterator(A2,n,kd)[0]]
    print("final amount2:",amount2)
    return bla,amount
'''
    #return

def getcyclic():
    while True:
        G=zero_matrix(F,0,n)
        aux=[]
        for i in range(n):
            aux+=[FL[randrange(q)]]
        for j in range(kd):
            G=G.stack(vector(aux))
            aux=shift(aux)
        #print("C:",G)
        if rank(G)==kd:
            #print("won")
            return G
        #print("lost")
        
    
        
def shift(lst):
    return [lst[-1]]+lst[:-1]

def getRM(r,m):
    #print("getting")
    if r==0:
        return matrix(F,[1]*(2**m))
    if r==m:
        return identity_matrix(F,2**m)
    G1=getRM(r,m-1)
    G2=getRM(r-1,m-1)
    G=G1.augment(G1)
    G2=zero_matrix(F,G2.nrows(),2**(m-1)).augment(G2)
    G=G.stack(G2)
    return G



mg=7
FF = GF(2**mg,'b')
b=FF.gen()
#PFF.<x> = PolynomialRing(FF) #python
def goppaiterator(t1):
    assert t1*mg<=n
    assert n<=2**mg
    #m
    #FF = GF(2^m,'a')
    #a=FF.gen()
    amount=[]
    bla=[]
    yaf24=0
    F1=FF.list()
    #D1=len(F1)
    D1=2**mg
    #D1=q?
    D2=t1+1
    S1=Subsets(range(D1),n)
    for l in range(1,D2):
        G1=Subsets(range(t1),l)
        for gp in G1:
            #print("gp:",gp)
            gp2=[0]*t1
            for v in gp:
                gp2[v]=1
            g=gp2+[1]
            g2=0
            for j in range(D2):
                g2+=(x**j)*g[j]
            #print("g2:",g2,g2.is_irreducible()) #if t=1 we can have zeroes; for now t>1 or it breaks
            #print("g:",g2,type(g2))
            if not g2.is_irreducible():
                continue
            for s in S1:
                L1=[]
                for z in s:
                    L1+=[F1[z]]
                A=mGoppa(g,L1).rref()
                while A[-1]==vector([0]*n): #dont know where to place
                    A=A[:-1]
                #print("GOPPA:",A)
                #print("L1:",L1)
                
                A=copy(A)
                global kd #can do like that?
                kd=rank(A)
                print("k:",kd)
                prip=SSAprep(dsf(A)[1])
                for h in bla:
                    A2=copy(A)#deepcopy
                    #prap=deepcopy(prip)
                    ess=SSA(h,prip,A2)
                    if type(ess)==list:
                        amount[bla.index(h)]+=1
                        yaf24=1
                        break
                if yaf24==1:
                    yaf24=0
                    continue
            #print("newclass")
                #print("adding A:",A)
                #A[0,0]=1
                print("adding")
                bla+=[copy(A)]#deepcopy
                amount+=[1]
    print("blaamount:",bla,amount)
    #print("nk",n,kd)
    amount2=[]
    for A in bla:
        A2=copy(A)#deepcopy
        kd=A2.nrows()
        amount2+=[iterator(A2,n,kd)[0]]
    print("final amount2:",amount2)
    return bla,amount

def goppastat(t1):
    assert t1*mg<=n
    assert n<=2**mg
    
    amount=[]
    bla=[]
    yaf24=0
    F1=FF.list()
    #D1=len(F1)
    D1=2**mg
    #D1=q?
    #D2=t1+1
    S1=Subsets(range(D1),n)
    G1=Subsets(range(t1),t1+1)
    hk=[0]*kd
    hfl=0
    for _ in range(1000000):
        if hfl==50:
            print("yay")
        if hfl==10000:
            break
        gp2=[0]*t1
        for v in range(t1):
            gp2[v]=F1[randrange(q)]
        g=gp2+[1]
        g2=0
        for j in range(t1+1):
            g2+=(x**j)*g[j]
        #print("g2:",g2,g2.is_irreducible()) #if t=1 we can have zeroes; for now t>1 or it breaks
        #print("g:",g2)
        if not g2.is_irreducible():
            #print("pluh")
            continue
        hfl+=1
        s=S1.random_element()
        L1=[]
        for z in s:
            L1+=[F1[z]]
        A=mGoppa(g,L1)
        #print(rank(A))
        if rank(A)!=mg*t1:
            print("pluh")
        hd=hully(dsf(A)[1])[0].nrows()
        hk[hd]+=1
        
    print(hk)

def LDPChelper(j,l):
    ar=[1]*l+[0]*(n-l)
    H=matrix(F,0,n)
    for i in range((n-kd)/j):
        H=H.stack(vector(ar))
        for k in range(l):
            ar=shift(ar)
    return H

def randLDPC(j,l):
    LD=matrix(F,0,n)
    H=LDPChelper(j,l)
    for i in range(j):
        perm=Permutations(n).random_element()
        pemr = Permutation(perm).to_matrix()
        HP=H*pemr
        LD=LD.stack(HP)
    return LD

#n,kd=7,3
#print(GRSiterator())
n,kd=128,64
#goppaiterator(6)
#goppastat(8)
#GRSiterator()

#getcyclic()
h0=0
h1=0
hk=[0]*kd
'''
for h in range(10000):
    if h==50:
        print("yay")
    C=getcyclic()
    #C=randLDPC(16,32)
    #C=random_matrix(F,kd,n)
    hd=hully(dsf(C)[1])[0].nrows()
    hk[hd]+=1
    #print(h0,h1)
    #if h0==10000:
     #   print(h0,h1)

print("hk:",hk)
'''

'''
m=10
C=getRM(5,m)
n=2**m
kd=C.nrows()
print("dim:",hully(dsf(C)[1])[0].nrows())
'''


def puncc2(C,pnc):
    #C6=copy(C)
    C6=C[:,:pnc].augment(C[:,pnc+1:])
    return C6
    
def dualc2(R):
    global n
    n=R.ncols()
    R=dsf(R)[1]
    l=-1
    H=zero_matrix(F,0,n)
    for k in range(n):
        if R[k,k]==0:
            
            H=H.stack(vector([0]*n))
            l+=1
            H[l,k]=1
            for j in range(k):
                H[l,j]=-R[j,k]
    return H
def shortc2(C,pnc):
    C=dsf(C)[1]
    Cp=copy(C)
    print("beforeshort:",pnc)
    print(Cp)
    l=Cp.nrows()
    for k in reversed(range(l)):
        if C[k,pnc]:
            print("k:",k)
            for i in range(k):
                if C[i,pnc]:
                    print("pang",i)
                    Cp[i]=C[i]+C[k]
            print("predelete")
            print(Cp)
            Cp=Cp.delete_rows([k])
            print("short before pnc:")
            print(Cp)
            #C6=Cp[:,:pnc].augment(Cp[:,pnc+1:])
            return Cp
    #C6=Cp[:,:pnc].augment(Cp[:,pnc+1:])
    return Cp

'''
Cpr = random_matrix(F,kd,n)
Cpr=dsf(Cpr)[1]
pnc=[randrange(n),randrange(n)]
print("Cpr:")
#print(Cpr)
#print(puncc2(Cpr,pnc))
print(Cpr)
print("pnc:")
print(pnc)
print("1:")
C2=dualc2(Cpr)
print("dualcpr:")
print(C2)
print("1.1")
for p in pnc:
    C2=shortc2(C2,p)
print("short:")
print(C2)
print("1.2")
C2=dsf(C2)[1]
for p in pnc:
    C2[p,p]=1
print("1.3")

print(C2)
print("2:")
C3=copy(Cpr)
for p in pnc:
    C3=puncc(C3,p)
print("pc3:")
print(C3)
print("2.1")
C3=dualc2(C3)
print("dc3:")
print(C3)
print("2.2")
C3=dsf(C3)[1]
print("2.3")

print(C3)
print("C2")
print(C2)
print("C3")
print(C3)
assert C2==C3
'''


n,kd=128,64
whatweight=0
dualflag=0
GMM=[]
rnd=100
altref=1


def importedSSA(C1,C2):
    C2=dsf(C2)[1]
    pre=SSAprep(C2)
    return SSA(C1,pre,C2)[0]
    
#example:
'''
Cpr = random_matrix(F,kd,n)
Cpr=dsf(Cpr)[1]
prip=SSAprep(Cpr)
C3=copy(Cpr)
debugperm=Permutations(n).random_element()
#debugperm=[12,11,10,9,8,7,1,2,3,4,5,6]
debugpemr = Permutation(debugperm).to_matrix()
debugpemr2= Permutation(debugperm).inverse().to_matrix()
C3=C3*debugpemr
debugperm2=[i-1 for i in debugperm]
#prem=SSA(C3,prip,Cpr)[0]
prem=importedSSA(C3,Cpr)
if type(prem) is not str:
    pmre=[i+1 for i in prem]
    pemr = Permutation(pmre).inverse().to_matrix()
    print("final:",dsf(C3*pemr)[1]==Cpr)
else:
    print("result:",prem)
'''
multilinecommentstringdeleter=1
