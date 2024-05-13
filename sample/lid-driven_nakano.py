import  numpy as np
from matplotlib import pyplot
plot_every=10


# Simulation parameters
Nx          = 20   # resolution x-dir
Ny          = 20    # resolution y-dir
rho0        = 1.0    # average density
tau         = 0.56    # collision timescale
Nt          = 1000  # number of timesteps

# Lattice speeds / weights
NL = 9
cxs = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1])
cys = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1])
weights = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36]) # sums to 1
uw = np.array([0.01, 0.0]) # wall velocity

# Initial Conditions - flow to the right with some perturbations
ux=np.zeros((Ny,Nx))
uy=np.zeros((Ny,Nx))

rho=np.zeros((Ny,Nx))
rho[:,:]=rho0
F=np.ones((Ny,Nx,NL))

def collision(F,rho, ux, uy):
    Feq=np.zeros(F.shape)
    for i,cx,cy,w in zip(range(NL),cxs,cys,weights):
      Feq[:,:,i]= rho*w*(1+3*(cx*ux+cy*uy)+9*(cx*ux+cy*uy)**2/2-3*(ux**2+uy**2)/2)
    F=F-(1/tau)*(F-Feq)
    return F

def streaming(F):
   for j in range(Ny):
      for i in range(Nx-1,0,-1):
         F[j,i,1]=F[j,i-1,1]
      for i in range(Nx-1):
         F[j,i,3]=F[j,i+1,3]
   for j in range(Ny-1):
      for i in range(Nx):
         F[j,i,2]=F[j+1,i,2]
      for i in range(Nx-1,0,-1):
         F[j,i,5]=F[j+1,i-1,5]
      for i in range(Nx-1):
         F[j,i,6]=F[j+1,i+1,6]
   for j in range(Ny-1,0,-1):
      for i in range(Nx):
         F[j,i,4]=F[j-1,i,4]
      for i in range(Nx-1):
         F[j,i,7]=F[j-1,i+1,7]
      for i in range(Nx-1,0,-1):
         F[j,i,8]=F[j-1,i-1,7]
   return F

def sfbound(F):
   #x=0, x=Nx-1, y=0 is stationary wall.
   #y=Ny-1 is moving wall.
   
   #West boundary
   F[:,0,1]=F[:,0,3]
   F[:,0,5]=F[:,0,7]
   F[:,0,8]=F[:,0,6]
   #East boundary
   F[:,Nx-1,3]=F[:,Nx-1,1]
   F[:,Nx-1,6]=F[:,Nx-1,8]
   F[:,Nx-1,7]=F[:,Nx-1,5]
   #South Boundary
   F[0,:,2]=F[0,:,4]
   F[0,:,5]=F[0,:,7]
   F[0,:,6]=F[0,:,8]
   #moving north boundary
   rhomb=np.ones(Nx)
   for i in range(1,Nx-1,1):
      rhomb[i]=F[0,i,0]+F[0,i,1]+F[0,i,3]+2*(F[0,i,4]+F[0,i,7]+F[0,i,8]) #try 265 instead of 478
   
   F[Ny-1,:,4]=F[Ny-1,:,2] - rhomb*(cxs[2]*uw[0] + cys[2]*uw[1]) / 6.0
   F[Ny-1,:,7]=F[Ny-1,:,5] - rhomb*(cxs[5]*uw[0] + cys[5]*uw[1]) / 6.0
   F[Ny-1,:,8]=F[Ny-1,:,6] - rhomb*(cxs[6]*uw[0] + cys[6]*uw[1]) / 6.0
   return F

def rhouv(F):
   rho=np.sum(F,2)
   for i in range(Nx):
      rho[i]=F[0,i,0]+F[0,i,1]+F[0,i,3]+2*(F[0,i,4]+F[0,i,7]+F[0,i,8])
   ux=np.sum(F*cxs,2)/rho
   uy=np.sum(F*cys,2)/rho
   return rho, ux, uy

def resultoutput(ux,uy,it):
    if(it%plot_every==0):
      print(it)
      pyplot.imshow(np.sqrt(ux**2+uy**2))
      pyplot.pause(.01)
      pyplot.cla()
    return


#main loop
for i in range(Nt):
   F=collision(F,rho,ux,uy)
   F=streaming(F)
   F=sfbound(F)
   rho,ux,uy=rhouv(F)
   resultoutput(ux,uy,i)
