import numpy as np
from matplotlib import pyplot
plot_every=10

# Simulation parameters
Nx          = 20   # resolution x-dir
Ny          = 20    # resolution y-dir
rho0        = 1.0    # average density
tau         = 0.56    # collision timescale
Nt          = 100000  # number of timesteps

X, Y = np.meshgrid(range(Nx), range(Ny))

# Lattice speeds / weights
NL = 9
cxs = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
cys = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
weights = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36]) # sums to 1


def initial_para():
   p=np.ones((Ny,Nx))
   ux=np.zeros((Ny,Nx))
   uy=np.zeros((Ny,Nx))
   return p, ux, uy

def bound_con(p,ux,uy):
   ux[Ny-1,1:Nx-1]=0.1
   uy[Ny-1,:]=0
   uy[0,:]=abs(uy[0,:])
   ux[:,0]=abs(ux[:,0])
   ux[:,Nx-1]=-abs(ux[:,Nx-1])
   #set boundary pressure to be same
   '''
   p_temp=1/3*(np.mean(p[0,:])+np.mean(p[:,0])+np.mean(p[:,Nx-1]))
   p[0,:]=p_temp
   p[1,:]=p_temp
   p[:,0]=p_temp
   p[:,1]=p_temp
   p[:,Nx-1]=p_temp
   p[:,Nx-2]=p_temp
   '''

   return p, ux, uy

def new_p(p,ux,uy):
   A_p=0.01
   for n in range(2):
      Feq_edit=np.zeros((Ny,Nx,NL))
      for i,cx,cy,w in zip(range(NL),cxs,cys,weights):
         for x in range (Nx):
            for y in range (Ny):
               if x==0:
                  Feq_edit[y,x,1]=Feq_edit[y,x,3]
                  Feq_edit[y,x,5]=Feq_edit[y,x,7]
                  Feq_edit[y,x,8]=Feq_edit[y,x,6]
               if x==Nx-1:
                  Feq_edit[y,x,3]=Feq_edit[y,x,1]
                  Feq_edit[y,x,7]=Feq_edit[y,x,5]
                  Feq_edit[y,x,6]=Feq_edit[y,x,8]
               if y==0:
                  Feq_edit[y,x,2]=Feq_edit[y,x,4]
                  Feq_edit[y,x,5]=Feq_edit[y,x,7]
                  Feq_edit[y,x,6]=Feq_edit[y,x,8]
               if y==Ny-1:
                  Feq_edit[y,x,4]=Feq_edit[y,x,2]
                  Feq_edit[y,x,7]=Feq_edit[y,x,5]-6*1/36*ux[y,x]
                  Feq_edit[y,x,8]=Feq_edit[y,x,6]+6*1/36*ux[y,x]
               else:
                  x_index=(x-cx)%Nx
                  y_index=(y-cy)%Ny
                  Feq_edit[y,x,i]=w*(3*p[y_index,x_index]+3*(cx*ux[y_index,x_index]+cy*uy[y_index,x_index])+9*(cx*ux[y_index,x_index]+cy*uy[y_index,x_index])**2/2-3*(ux[y_index,x_index]**2+uy[y_index,x_index]**2)/2)
                  +3*A_p*w*(cx*(ux[y,x]-ux[y_index,x_index])+cy*(uy[y,x]-uy[y_index,x_index]))
      
      p=np.sum(Feq_edit,2)/3
   return p

def new_u(p,ux,uy):
   A_u=0.01
   Feq_edit=np.zeros((Ny,Nx,NL))
   for i,cx,cy,w in zip(range(NL),cxs,cys,weights):
      for x in range (Nx):
         for y in range (Ny):
               if x==0:
                  Feq_edit[y,x,1]=Feq_edit[y,x,3]
                  Feq_edit[y,x,5]=Feq_edit[y,x,7]
                  Feq_edit[y,x,8]=Feq_edit[y,x,6]
               if x==Nx-1:
                  Feq_edit[y,x,3]=Feq_edit[y,x,1]
                  Feq_edit[y,x,7]=Feq_edit[y,x,5]
                  Feq_edit[y,x,6]=Feq_edit[y,x,8]
               if y==0:
                  Feq_edit[y,x,2]=Feq_edit[y,x,4]
                  Feq_edit[y,x,5]=Feq_edit[y,x,7]
                  Feq_edit[y,x,6]=Feq_edit[y,x,8]
               if y==Ny-1:
                  Feq_edit[y,x,4]=Feq_edit[y,x,2]
                  Feq_edit[y,x,7]=Feq_edit[y,x,5]-6*1/36*ux[y,x]
                  Feq_edit[y,x,8]=Feq_edit[y,x,6]+6*1/36*ux[y,x]
               else:
                  x_index=(x-cx)%Nx
                  y_index=(y-cy)%Ny
                  Feq_edit[y,x,i]=w*(3*p[y_index,x_index]+3*(cx*ux[y_index,x_index]+cy*uy[y_index,x_index])+9*(cx*ux[y_index,x_index]+cy*uy[y_index,x_index])**2/2-3*(ux[y_index,x_index]**2+uy[y_index,x_index]**2)/2)
                  +3*A_u*w*(cx*(ux[y,x]-ux[y_index,x_index])+cy*(uy[y,x]-uy[y_index,x_index]))
   ux=np.sum(Feq_edit*cxs,2)
   uy=np.sum(Feq_edit*cys,2)
   
   return ux,uy

def resultoutput(ux,uy,it):
    if(it%plot_every==0):
      print(it)
      
      pyplot.quiver(X,Y,ux*20,uy*20,color=(1, 0, 0, 0.3),angles='xy', scale_units='xy', scale=1)
      pyplot.grid()
      pyplot.show()
      '''
      pyplot.imshow(np.sqrt(ux**2+uy**2))
      pyplot.pause(.01)
      pyplot.cla()
      '''
    return



#main loop
p,ux, uy =initial_para()
p,ux,uy=bound_con(p,ux,uy)
for t in range(Nt):
   #p,ux,uy=bound_con(p,ux,uy)
   p=new_p(p,ux,uy)
   ux,uy=new_u(p,ux,uy)
   p,ux,uy=bound_con(p,ux,uy)
   resultoutput(ux,uy,t)