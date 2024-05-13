import numpy as np
from matplotlib import pyplot
import math
plot_every=10
# Simulation parameters
Nx          = 40   # resolution x-dir
Ny          = 40    # resolution y-dir
rho0        = 1.0    # average density
Nt          = 1600  # number of timesteps

X, Y = np.meshgrid(range(Nx), range(Ny))

# Lattice speeds / weights
NL = 9
cxs = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
cys = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
weights = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36]) # sums to 1
A = 0.0
uw = np.array([0.1, 0])
k1 = 2*3.14/Nx
k2 = 2*3.14/Ny
nu = 1/6 - (2*A)/9


def initial_para():
    p=np.zeros((Ny,Nx))
    ux=np.zeros((Ny,Nx))
    uy=np.zeros((Ny,Nx))
    feq=np.zeros((Ny,Nx,NL))
    f=np.zeros((Ny,Nx,NL))
    grad_u=np.zeros((Ny,Nx,2,2))
    cr=np.zeros((NL,2))
    
    for x in range (Nx):
        for y in range (Ny):
            ux[y,x] = -0.01*math.cos(k1*x)*math.sin(k2*y)
            uy[x,y] = k1/k2*0.01*math.sin(k1*x)*math.cos(k2*y)
            p[:,:]= 1/3 - 0.01*0.01*0.25*(math.cos(2*k1*x)+(k1/k2)**2*math.cos(2*k2*y))
    return p, ux, uy, feq, f, grad_u, cr

def cr_cal(cr):
    for i in range (NL):
        cr[i,0]=cxs[i]
        cr[i,1]=cys[i]
    
    return cr

def bound_con(p,ux,uy):
    ux[Ny-1,:]=uw[0]
    uy[Ny-1,:]=0
    uy[0,:]=0
    ux[:,0]=0
    ux[:,Nx-1]=0

    p[:,0]=p[:,1]
    p[:,Nx-1]=p[:,Nx-2]
    p[0,:]=p[1,:]

    return p, ux, uy

def grad_u_cal(grad_u,ux,uy,cr):
    for x in range (Nx):
        for y in range (Ny):
            for beta in range(2):
                grad_u[y,x,beta,0]=0
                grad_u[y,x,beta,1]=0
                for i in range (NL):
                    x_index=(x+cxs[i])%Nx
                    y_index=(y+cys[i])%Ny
                    grad_u[y,x,beta,0]=grad_u[y,x,beta,0] + cr[i,beta]*ux[y_index,x_index]
                    grad_u[y,x,beta,1]=grad_u[y,x,beta,1] + cr[i,beta]*uy[y_index,x_index]
                grad_u[y,x,beta,0]=grad_u[y,x,beta,0]/6
                grad_u[y,x,beta,1]=grad_u[y,x,beta,1]/6
                
    return grad_u

def feq_cal(feq,ux,uy,p,cr):
    for i in range (NL):
        for x in range (Nx):
            for y in range (Ny):
                ftemp = 0
                for alpha in range (2):
                    for beta in range (2):
                        ftemp = ftemp + (grad_u[y,x,beta,alpha]+grad_u[y,x,alpha,beta])*cr[i,alpha]*cr[i,beta]
                feq[y,x,i]=weights[i]*(3*p[y,x]+3*(cxs[i]*ux[y,x]+cys[i]*uy[y,x])+4.5*(cxs[i]*ux[y,x]+cys[i]*uy[y,x])**2-1.5*(ux[y,x]*ux[y,x]+uy[y,x]*uy[y,x]))
                
                feq[y,x,i]=feq[y,x,i] + A*weights[i]*ftemp
                
    return feq

def physics(p,ux,uy,feq):
    ux[:,:]=0
    uy[:,:]=0
    p[:,:]=0
    for x in range (Nx):
        for y in range (Ny):
            for i in range (NL):
                x_index=(x-cxs[i])%Nx
                y_index=(y-cys[i])%Ny
                ux[y,x]=ux[y,x] + cxs[i]*feq[y_index,x_index,i]
                uy[y,x]=uy[y,x] + cys[i]*feq[y_index,x_index,i]
                p[y,x]=p[y,x] + feq[y_index,x_index,i]/3
                
    return p,ux,uy

def bounceback(feq):
    y = 0
    for x in range (Nx):
        x_index=(x-cxs[2])%Nx
        y_index=(y-cys[2])%Ny
        feq[y_index,x_index,2]=feq[y_index,x_index,4]
        x_index=(x-cxs[5])%Nx
        y_index=(y-cys[5])%Ny
        feq[y_index,x_index,5]=feq[y_index,x_index,7]
        x_index=(x-cxs[6])%Nx
        y_index=(y-cys[6])%Ny
        feq[y_index,x_index,6]=feq[y_index,x_index,8]
        
    y = Ny-1
    for x in range (Nx):
        x_index=(x-cxs[4])%Nx
        y_index=(y-cys[4])%Ny
        feq[y_index,x_index,4]=feq[y_index,x_index,2] - 6*weights[2]*(cxs[2]*uw[0]+cys[2]*uw[1])
        x_index=(x-cxs[7])%Nx
        y_index=(y-cys[7])%Ny
        feq[y_index,x_index,7]=feq[y_index,x_index,5] - 6*weights[5]*(cxs[5]*uw[0]+cys[5]*uw[1])
        x_index=(x-cxs[8])%Nx
        y_index=(y-cys[8])%Ny
        feq[y_index,x_index,8]=feq[y_index,x_index,6] - 6*weights[6]*(cxs[6]*uw[0]+cys[6]*uw[1])
        
    x = 0
    for y in range (Ny):
        x_index=(x-cxs[1])%Nx
        y_index=(y-cys[1])%Ny
        feq[y_index,x_index,1]=feq[y_index,x_index,3]
        x_index=(x-cxs[5])%Nx
        y_index=(y-cys[5])%Ny
        feq[y_index,x_index,5]=feq[y_index,x_index,7]
        x_index=(x-cxs[8])%Nx
        y_index=(y-cys[8])%Ny
        feq[y_index,x_index,8]=feq[y_index,x_index,6]
    
    x = Nx-1
    for y in range (Ny):
        x_index=(x-cxs[3])%Nx
        y_index=(y-cys[3])%Ny
        feq[y_index,x_index,3]=feq[y_index,x_index,1]
        x_index=(x-cxs[7])%Nx
        y_index=(y-cys[7])%Ny
        feq[y_index,x_index,7]=feq[y_index,x_index,5]
        x_index=(x-cxs[6])%Nx
        y_index=(y-cys[6])%Ny
        feq[y_index,x_index,6]=feq[y_index,x_index,8]
        
    return feq

def bounceback2(feq):
    y = 0
    for x in range (Nx):
        feq[y,x,2]=feq[y,x,4]
        feq[y,x,5]=feq[y,x,7]
        feq[y,x,6]=feq[y,x,8]
        
    y = Ny-1
    for x in range (Nx):
        feq[y,x,4]=feq[y,x,2] - 6*weights[2]*(cxs[2]*uw[0]+cys[2]*uw[1])
        feq[y,x,7]=feq[y,x,5] - 6*weights[5]*(cxs[5]*uw[0]+cys[5]*uw[1])
        feq[y,x,8]=feq[y,x,6] - 6*weights[6]*(cxs[6]*uw[0]+cys[6]*uw[1])
        
    x = 0
    for y in range (Ny):
        feq[y,x,1]=feq[y,x,3]
        feq[y,x,5]=feq[y,x,7]
        feq[y,x,8]=feq[y,x,6]
    
    x = Nx-1
    for y in range (Ny):
        feq[y,x,3]=feq[y,x,1]
        feq[y,x,7]=feq[y,x,5]
        feq[y,x,6]=feq[y,x,8]
        
    return feq

def physics2(p,ux,uy,feq):
    ux[:,:]=0
    uy[:,:]=0
    p[:,:]=0
    for x in range (Nx):
        for y in range (Ny):
            for i in range (NL):
                x_index=(x+cxs[i])%Nx
                y_index=(y+cys[i])%Ny
                ux[y_index,x_index]=ux[y_index,x_index] + cxs[i]*feq[y,x,i]
                uy[y_index,x_index]=uy[y_index,x_index] + cys[i]*feq[y,x,i]
                p[y_index,x_index]=p[y_index,x_index] + feq[y,x,i]/3
                
    return p,ux,uy

def resultoutput(ux,uy,it):
    if(it%plot_every==0):
      print(it)
      
      pyplot.quiver(X,Y,ux*100,uy*100,color=(1, 0, 0, 0.3),angles='xy', scale_units='xy', scale=1)
      pyplot.grid()
      pyplot.show()
      '''
      pyplot.imshow(np.sqrt(ux**2+uy**2))
      pyplot.pause(.01)
      pyplot.cla()
      '''
    return



#main loop
p,ux,uy,feq,f,grad_u,cr =initial_para()
cr = cr_cal(cr)

for t in range(Nt):
    
    grad_u=grad_u_cal(grad_u,ux,uy,cr)
    feq=feq_cal(feq,ux,uy,p,cr)
    # feq=bounceback(feq)
    p,ux,uy=physics(p,ux,uy,feq)
    
    resultoutput(ux,uy,t)

#    p=new_p(p,ux,uy)
#    ux,uy=new_u(p,ux,uy)
#    p,ux,uy=bound_con(p,ux,uy)
    