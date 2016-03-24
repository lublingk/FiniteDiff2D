# Define parameters:
xleft = 0
xright = 2*pi
yleft = 0
yright = pi
a = 1
b = 1
mx = 50
my = 54
def f(x,y):
    return 2*a*sin(x)*sin(y)+b*(sin(x)*sin(y)+x+y)
  
def g(x,y):
    return x+y

exact_sol_defined = True

def exact(x,y):
    return sin(x)*sin(y)+x+y

fdm2d(xleft, xright, yleft, yright, a, b, g, g, g, g, f, mx, my, exact_sol_defined, exact)

from numpy import *
from scipy.linalg import solve
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

import pylab as plt

# xleft, xright ... interval
# a, b ... constant coefficients
# ua, ub ... Dirichlet conditions
# f ... right-hand side
# m ... number of subdivisions
def fdm2d(xleft, xright, yleft, yright, a, b, uxl, uxr, uyl, uyr, f, mx, my, exact_sol_defined = False, exact = []):
    # Calculate h and k:
    h = (xright - xleft) / float(mx)
    k = (yright - yleft) / float(my)
    # Size of rhs vector:
    n = (mx - 1)*(my - 1)
    # Predefine empty coordinate matrix:
    M = []
    Mrow = []
    Mcol = []
    # Right-hand side vector:
    F = zeros(n)
    
    #Calculate matrix
    
    #initialize count
    count = 0
    
    #Loop through all interior points
    for I in range(mx-1):
      for J in range(my-1):

        F[count]=(f(xleft+(I+1)*h,yleft+(J+1)*k))
        
        #If primary x left interior point
        if(I == 0):
          F[count]=F[count] + (a*uxl(xleft+I*h,yleft+(J+1)*k))/(h**2)
          
        else:
          M.append((-a)/(h**2))
          Mrow.append(count)
          Mcol.append(count-my+1)

        #If primary y left interior point
        if(J == 0):
          F[count]=F[count] + (a*uyl(xleft+(I+1)*h,yleft+J*k))/(k**2)
          
        else:
          M.append((-a)/(k**2))
          Mrow.append(count)
          Mcol.append(count-1)
                
        M.append((2*a/(h**2))+(2*a/(k**2))+b)
        Mrow.append(count)
        Mcol.append(count)
          
        #If primary y right interior point
        if(J == (my - 2)):
          F[count]=F[count] + (a*uyr(xleft+(I+1)*h,yleft+(J+2)*k))/(k**2)

        else:
          M.append((-a)/(k**2))
          Mrow.append(count)
          Mcol.append(count+1)
        
        
        #If primary x right interior point
        if(I == (mx - 2)):
          F[count]=F[count] + (a*uxr(xleft+(I+2)*h,yleft+(J+1)*k))/(h**2)
          
        else:
          M.append((-a)/(h**2))
          Mrow.append(count)
          Mcol.append(count+my-1)
          
        count=count + 1
        
    Mcsr=csr_matrix( (M,(Mrow,Mcol)), shape=(n,n) ).tocsr()
    
    resultvect=spsolve(Mcsr,F)
    pdmresult=zeros((mx-1,my-1))
    count=0
    for I in range(mx-1):
      for J in range(my-1):
        pdmresult[I][J]=resultvect[count]
        
        count=count+1

    coordx=zeros(mx+1)
    coordy=zeros(my+1)
    plotresult=zeros((mx+1,my+1))
    
    if exact_sol_defined:
      exactcoordx=zeros(2*mx+1)
      exactcoordy=zeros(2*my+1)
      exactresult=zeros((2*mx+1,2*my+1))
      for I in range(2*mx+1):
        exactcoordx[I]=xleft+I*float(h/2)
      for J in range(2*my+1):
        exactcoordy[J]=yleft+J*float(k/2)
        
      for I in range(2*mx+1):
        for J in range(2*my+1):
          exactresult[I][J]=exact(xleft+I*float(h/2),yleft+J*float(k/2))
          
      

    for I in range(mx+1):    
      coordx[I]=xleft+I*h
      
    for J in range(my+1):
      coordy[J]=yleft+J*k
      
    for I in range(mx+1):
      for J in range(my+1):
        if(I==0):
          plotresult[I][J]=uxl(xleft+I*h,yleft+J*k)
          
        elif(J==0):
          plotresult[I][J]=uyl(xleft+I*h,yleft+J*k)
          
        elif(I==mx):
          plotresult[I][J]=uxr(xleft+I*h,yleft+J*k)
          
        elif(J==my):
          plotresult[I][J]=uyr(xleft+I*h,yleft+J*k)
          
        else:
          plotresult[I][J]=pdmresult[I-1][J-1]


    coordx,coordy = meshgrid(coordy,coordx)
          
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.plot_surface(coordy,coordx,plotresult,rstride=2,cstride=2)
    ax.set_xlabel('x-axis')
    ax.set_ylabel('y-axis')
    ax.set_zlabel('z-axis')
    
    lab.show()
    
    if exact_sol_defined:
      
      exactcoordx,exactcoordy=meshgrid(exactcoordy,exactcoordx)
      
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')

      surf = ax.plot_surface(exactcoordy,exactcoordx,exactresult,rstride=2,cstride=2)
      ax.set_xlabel('x-axis')
      ax.set_ylabel('y-axis')
      ax.set_zlabel('z-axis')
    
      lab.show()
    