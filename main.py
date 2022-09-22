#reference https://ieeexplore-ieee-org.sproxy.dongguk.edu/stamp/stamp.jsp?tp=&arnumber=4357612
import sympy as sym
import math
import numpy as np
import numpy.linalg as LA
import random
#Pk -> tk 초기화
#try3
Pk=[[0 for k in range(4)] for j in range(3)]#power spectra
tk=[[0 for k in range(4)] for j in range(3)] # user k가 다른 사용자에게 미치는 간섭 영향
alpha=[1 for k in range(3)]


K=3
N=4
Ik={0:[0,0,0,0],1:[0,0,0,0],2:[0,0,0,0]}

d = [[3, 20, 21],
     [20, 4, 22],
     [21, 22, 5]]

H=dict()
for i in range(3):
    p=[[0 for k in range(4)] for j in range(3)]
    H[i]=p
    

for i in range(K):
    for j in range(K):
        for k in range(N):
            if i == j:  #데이터 채널
                H[i][j][k] = d[i][i]**(-2) * random.expovariate(1)  #path loss * fading
            else:   #간섭 채널
                H[i][j][k] = d[i][j] ** (-2) * random.expovariate(1)   #path loss * fading


noise=1e-7


for key in range(3):
    p=np.random.randint(low=0,high=10,size=4)
    p=np.array(p/np.sum(p))*100
    Pk[key]=p





#tk(n) 식
for i in range(K):
    for j in range(N): #index
        tmp_0 = 0
        for k in range(K): 
            tmp = 0
            for l in range(K):
                if l != k:
                    tmp = tmp + Pk[l][j] * H[l][k][j]
            if k != i:
                tmp_0 = tmp_0 + (alpha[k]*Pk[k][j]*H[k][k][j] / (Pk[k][j] * H[k][k][j] + tmp + noise)) * (H[i][k][j] / (tmp + noise))
        tk[i][j] = tmp_0
#k값 default

    

pmax=100
threshold=1e-3
lamb=[0 for k in range(3)]
lmin=[0 for k in range(3)]
lmax=[10 for k in range(3)]


def f(x):
    return pmax-x
loop1=1
loop2=1
loop3=1
while(loop1):
    while(loop2):
        for k in range(K):
            for n in range(N):
                temp_Ik = 0
                for j in range(K):
                    if j != k:
                        temp_Ik = temp_Ik + (Pk[j][n] * H[j][k][n])
                Ik[k][n] = temp_Ik
        
        px=[0 for i in range(3)] #lmin일 때 power값 합
        py=[0 for i in range(3)] #lmax일때 power값 합
        pnew=[0 for i in range(3)] # 새 람다로 구했을 때 power값합
        
        for k in range(3):#람다구하기
            lmin[k]=0
            lmax[k]=5
            loop3=1
            while(loop3): 
                #print(k,lamb[k])
                px[k]=0
                py[k]=0
                pnew[k]=0
                for n in range(4):
                    temp1=alpha[k]/(lmin[k]+tk[k][n])-(Ik[k][n]+noise)/H[k][k][n]
                    temp2=alpha[k]/(lmax[k]+tk[k][n])-(Ik[k][n]+noise)/H[k][k][n]
                    if temp1>0:
                        px[k]=px[k]+temp1
                    if temp2>0:
                        py[k]+=temp2
                #print('lmin,lmax p합',k,px[k],py[k])      
                if f(px[k])*f(py[k])>0:
                    loop3=0 # 해가 존재하지 않으면
                elif f(px[k])*f(py[k])<0:
                    lamb[k]=(lmax[k]+lmin[k])/2 # 이분법
                    for n in range(4):
                        temp3=alpha[k]/(lamb[k]+tk[k][n])-(Ik[k][n]+noise)/H[k][k][n]
                        if temp3>0:
                            pnew[k]+=temp3
                            
                    if f(pnew[k])<threshold and f(pnew[k])>0:
                        loop3=0 #해가 존재하면
                    elif f(pnew[k])>0:
                        lmax[k]=lamb[k]
                    elif f(pnew[k])<0:
                        lmin[k]=lamb[k]
                        
                        
        Pnew=[[0 for k in range(4)] for j in range(3)]
        for k in range(3):
            for n in range(4):
                temp=alpha[k]/(lamb[k]+tk[k][n])-(Ik[k][n]+noise)/H[k][k][n]
                if temp>0:
                    Pnew[k][n]=temp
                #else:
                    #Pnew[k][n]=0
        #pk 수렴 check
        conv_pk = [[0 for j in range(N)] for i in range(K)]

        for k in range(K):
            for n in range(N):
                if math.sqrt((Pnew[k][n] - Pk[k][n]) * (Pnew[k][n] - Pk[k][n])) <= 0.001:
                    conv_pk[k][n] = 1

        conv_pk_check = 1
        for k in range(K):
            for n in range(N):
                if conv_pk[k][n] == 0:
                    conv_pk_check = 0
        for k in range(K):
            Pk[k]=Pnew[k]
        if conv_pk_check == 1:
            loop2 = 0
        #for k in range(K):
            #print('user',k,Pk[k])
                    

            
        
    tnew=[[0 for k in range(4)] for j in range(3)]
    for i in range(K):
        for j in range(N): #index
            tmp_0 = 0
            for k in range(K): 
                tmp = 0
                for l in range(K):
                    if l != k:
                        tmp = tmp + Pk[l][j] * H[l][k][j]
                if k != i:
                        tmp_0 = tmp_0 + (alpha[k]*Pk[k][j]*H[k][k][j] / (Pk[k][j] * H[k][k][j] + tmp + noise)) * (H[i][k][j] / (tmp + noise))
            tnew[i][j] = tmp_0
   #tk(n) 수렴 check
    conv_tk = [[0 for j in range(N)] for i in range(K)]

    for k in range(K):
        for n in range(N):
            if math.sqrt((tnew[k][n] - tk[k][n]) * (tnew[k][n] - tk[k][n])) <= 0.001:
                #print('check')
                conv_tk[k][n] = 1

    conv_tk_check = 1
    for k in range(K):
        for n in range(N):
            if conv_tk[k][n] == 0:
                conv_tk_check = 0
    for k in range(K):
        tk[k]=tnew[k]
    if conv_tk_check == 1:
        #print('tk빠짐')
        loop1 = 0
    

                    
print('mwif result')
for k in range(3):
    print('user',k,Pk[k],'pk합',np.sum(Pk[k]))
   # print(k,tk[k])

        
print('람다값:',lamb)

#capacity
capacity1=0
for k in range(3):
    for n in range(4):
        num=0
        for j in range(3):
            if j==k:
                continue
            num+=Pk[j][n]*H[j][k][n]
        capacity1=capacity1+alpha[k]*math.log10(1+(Pk[k][n]*H[k][k][n])/(num+noise))
        
print(capacity1)

#verifying

max=0
Preal=[[0 for k in range(N)] for j in range(K)]
for x in range(1000000):
    #ptemp 각각 합이 100인
    ptemp=[[0 for k in range(N)] for j in range(K)]
    for k in range(K):
        temp = 0
        for n in range(N):
            ptemp[k][n] = random.random()
            temp = temp + ptemp[k][n]

        for a in range(N):
            ptemp[k][a] = ptemp[k][a] * (100 / temp)

    capacity2=0
    for k in range(3):
        for n in range(4):
            num=0
            for j in range(3):
                if j==k:
                    continue
                num+=ptemp[j][n]*H[j][k][n]
            capacity2+=alpha[k]*math.log10(1+(ptemp[k][n]*H[k][k][n])/(num+noise))
    if capacity2>max:
        max=capacity2
        for i in range(3):
            Preal[i]=ptemp[i]
print('최적해')
for i in range(3):
    print('preal',i,Preal[i])
print(max)