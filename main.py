#!/usr/bin/python3
import numpy as np
from numpy import tile
from scipy import linalg as LA
import math,cmath
pi = math.pi
cos = math.cos
sin = math.sin
exp=np.exp
sqrt = math.sqrt
log=np.log



#----------------------------------------------------------------------------------------------------------
def band(t1,t2,a,Nk,density,mass,dim,totomega,m,sigma,square=True):
    print (t1,t2,a,Nk)
    if square:
        assert (t1 > 0), "Only positive hopping are possible!"
        print('Band structure plot of the square lattice with nearest and next nearest neighbour hopping bl......')
        gap1=np.zeros((Nk,Nk))
        gap2=np.zeros((Nk,Nk))
        gap3=np.zeros((Nk,Nk))
        gap4=np.zeros((Nk,Nk)) 
        basis=4
        wavefn_0=np.zeros((((Nk,Nk,basis,basis))))
        wavefn_1=np.zeros((((Nk,Nk,basis,basis))))
        wavefn_2=np.zeros((((Nk,Nk,basis,basis))))
        wavefn_3=np.zeros((((Nk,Nk,basis,basis))))
        f = open("square.dat", 'w')
        g = open("square_gamma_x_m_gamma.dat", 'w')
        h=  open('wavefn.dat','w');
        V=mass
        n=density
        alpha=0.0
        kx_m=[]
        for nkx in range(0, Nk):
            kx = (2.0*pi)/(Nk)*nkx
            kx_m.append(kx)
            ky_m=[]
            for nky in range(0, Nk):
                ky = (2.0*pi)/(Nk)*nky   
                ky_m.append(ky)
                #fk = 1.0+exp(-1j*kx*a)
                fk=cos(0.5*a*kx)
                #gk = 1.0+exp(-1j*ky*a)
                gk=cos(0.5*a*ky)
                
                #Ek = np.array([[0,(1+alpha)*t1*fk,(1+alpha)*t1*gk,0],[(1+alpha)*t1*np.conjugate(fk),0.0,0,(1-alpha)*V],[(1+alpha)*t1*np.conjugate(gk),0,0.0,(1-alpha)*V],[0,(1-alpha)*V,(1-alpha)*V,0.0]])
#                Ek = np.array([[0.0,(1+alpha)*t1*fk,-(1+alpha)*t1*gk,0],[(1+alpha)*t1*np.conjugate(fk),0.0,0,(1-alpha)*t1*gk],[-(1+alpha)*t1*np.conjugate(gk),0,0.0,(1-alpha)*t1*fk],[0,(1-alpha)*t1*np.conjugate(gk),(1-alpha)*t1*np.conjugate(fk),0.0]])
                Ek = np.array([[0,(1+alpha)*t1*fk,(1+alpha)*t1*gk,0],[(1+alpha)*t1*np.conjugate(fk),0.0,0,1.0*(1-alpha)*t1*gk],[(1+alpha)*t1*np.conjugate(gk),0,0.0,1.0*(1-alpha)*t1*fk],[0,1.0*(1-alpha)*t1*np.conjugate(gk),1.0*(1-alpha)*t1*np.conjugate(fk),-sigma*m*V]])
                #vk=LA.eigvals(Ek)
                vk, wk=np.linalg.eigh(Ek)
                #vk=sorted(vk)
                gap1[nkx][nky]=vk[0].real
                gap2[nkx][nky]=vk[1].real
                gap3[nkx][nky]=vk[2].real
                gap4[nkx][nky]=vk[3].real
                for rnbasis in range(0, basis):
                    for cnbasis in range(0, basis):
                        wavefn_0[nkx][nky][rnbasis][cnbasis]= np.dot(wk[rnbasis,0],wk[cnbasis,0]).real
                        wavefn_1[nkx][nky][rnbasis][cnbasis]= np.dot(wk[rnbasis,1],wk[cnbasis,1]).real
                        wavefn_2[nkx][nky][rnbasis][cnbasis]= np.dot(wk[rnbasis,2],wk[cnbasis,2]).real
                        wavefn_3[nkx][nky][rnbasis][cnbasis]= np.dot(wk[rnbasis,3],wk[cnbasis,3]).real
                        h.write("%5.2f %5.2f %5.2f %5.2f\n" % (wavefn_0[nkx][nky][rnbasis][cnbasis],wavefn_1[nkx][nky][rnbasis][cnbasis],wavefn_2[nkx][nky][rnbasis][cnbasis],wavefn_3[nkx][nky][rnbasis][cnbasis]))
                f.write("%5.2f %5.2f %5.2f %5.2f\n" % (gap1[nkx][nky],gap2[nkx][nky],gap3[nkx][nky],gap4[nkx][nky]))
                if((nky==0) & (nkx <= 0.5*Nk)):
                    g.write("%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (kx, ky, kx+ky, vk[0].real, vk[1].real, vk[2].real, vk[3].real))
        f.close()
        g.close()
        h.close()
       
        np.savetxt('kvecs.dat',np.c_[kx_m,ky_m],fmt='%1.4f %1.4f',header="# kx, ky",comments='')
        
        g = open("square_gamma_x_m_gamma.dat", 'a')

        for nkx in range(0, Nk):
            kx = (2*pi)/(Nk)*nkx
            kx_m.append(kx)
            ky_m=[]
            for nky in range(0, Nk):
                ky = (2*pi)/(Nk)*nky   

                if((nkx==(int(0.5*Nk)+1)) & (nky <= (int(0.5*Nk)))):
    			
                    g.write("%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (kx, ky, kx+ky,gap1[nkx][nky],gap2[nkx][nky],gap3[nkx][nky],gap4[nkx][nky]))
      			
        g.close()
        s = open("square_gamma_x_m_gamma.dat", 'a')



        for nkx in range((int(0.5*Nk)+1), Nk):
            kx = (2*pi)/(Nk)*nkx
            kx_m.append(kx)
            ky_m=[]
            for nky in range((int(0.5*Nk)+1), Nk):
                ky = (2*pi)/(Nk)*nky   

                if(nkx==nky):
                    s.write("%5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (kx, ky, kx+ky,gap1[nkx][nky],gap2[nkx][nky],gap3[nkx][nky],gap4[nkx][nky]))
#---------------------DOS--------------------------------------------------------------------------------------------------------------------

        
        
        f0 = open('dos0.dat','w')
        f1 = open('dos1.dat','w')
        f2 = open('dos2.dat','w')
        f3 = open('dos3.dat','w')
        
        jrkky=np.zeros((basis,basis))
        dosbasis=np.zeros(((basis,basis,totomega)));
        for rnbasis in range(0,basis):
            for cnbasis in range(0,basis):
                jrkky[rnbasis,cnbasis]=0.0
                for Nomega in range(1,int(totomega)):
                    omega=10.0/totomega*(Nomega-totomega/2)
                    dossum1 =0.0
                    dossum2 =0.0
                    dossum3 =0.0
                    dossum4 =0.0

                    for nkx in range (0,Nk):
                        for nky in range(0,Nk):
                            z=omega+0.04j
                            A1=-(1.0/pi)*((wavefn_0[nkx][nky][rnbasis][cnbasis]/(z-gap1[nkx][nky]))).imag
                    
                            dossum1=dossum1+A1
                            A2=-(1.0/pi)*((wavefn_1[nkx][nky][rnbasis][cnbasis]/(z-gap2[nkx][nky]))).imag
                            
                            dossum2=dossum2+A2
                            
                            A3=-(1.0/pi)*((wavefn_2[nkx][nky][rnbasis][cnbasis]/(z-gap3[nkx][nky]))).imag
                           
                            dossum3=dossum3+A3
                            
                            A4=-(1.0/pi)*((wavefn_3[nkx][nky][rnbasis][cnbasis]/(z-gap4[nkx][nky]))).imag
                            
                            dossum4=dossum4+A4
                        
                    dosbasis[rnbasis,cnbasis,Nomega]=dossum1+dossum2+dossum3+dossum4
                    if ((rnbasis==0) and (cnbasis==0)):
                        f0.write("%6.2f %6.2f\n" %(omega, dosbasis[rnbasis,cnbasis,Nomega]/(Nk*Nk)))
                    elif((rnbasis==1) and (cnbasis==1)):
                        f1.write("%6.2f %6.2f\n" %(omega, dosbasis[rnbasis,cnbasis,Nomega]/(Nk*Nk)))
                    elif((rnbasis==2) and (cnbasis==2)):
                        f2.write("%6.2f %6.2f\n" %(omega, dosbasis[rnbasis,cnbasis,Nomega]/(Nk*Nk)))
                    elif((rnbasis==3) and (cnbasis==3)):
                        f3.write("%6.2f %6.2f\n" %(omega, dosbasis[rnbasis,cnbasis,Nomega]/(Nk*Nk)))
                   
           
        #-----------------------------------------------------------norm of the dos-------------------------------------------------------------------
      
        
        #----------------------------------------------------------------------------------------------------------------------------------------------	
        
        #------------------------------------------------------------------------------------------------------------------------------------------------
    
    else:
        print('Band structure plot of the Lieb lattice (a decorated lattice) with nearest and next nearest neighbour hopping......')

        assert ( t1 > 0), 'Only positive hoppin are allowed!'	
        gapplus=np.zeros((Nk,Nk))
        gapflat=np.zeros((Nk,Nk))
        gapminus=np.zeros((Nk,Nk))
        basis=3
        wavefn_0=np.zeros(((Nk,Nk,basis)))
        wavefn_1=np.zeros(((Nk,Nk,basis)))
        wavefn_2=np.zeros(((Nk,Nk,basis)))

        f = open('dispersion_lieb.dat','w');
        g=  open('wavefn.dat','w');
        kx_m=[]
        for nkx in range(0, Nk):
            kx = 2*pi/Nk*nkx
            kx_m.append(kx)
            ky_m=[]
            for nky in range(0, Nk):
                ky = 2*pi/Nk*nky
                ky_m.append(ky)
                fk = cos(0.5*a*kx)
                gk = cos(0.5*a*ky)
                fkk = 2*cos(0.5*a*kx)*cos(0.5*a*ky)
                gkk = 2*cos(0.5*a*kx)*cos(0.5*a*ky)
                fdim = 2*(1j)*sin(0.5*a*kx)
                gdim = 2*(1j)*sin(0.5*a*ky)
                    #Ek = -2.0*t1*np.array([[0.0, fk, gk],[fk, 0.0, 0.0],[gk, 0.0, 0.0]])-4.0*t2*np.array([[0.0, 0.0, 0.0],[0.0, 0.0, fkk],[0.0, gkk, 0.0]])+dim*np.array([[0.0, -fdim, -gdim],[fdim, 0.0, 0.0],[ gdim, 0.0, 0.0]])

                    #Ek = -2.0*t1*np.array([[0.0, fk, 0.0],[fk, 0.0, 0.0],[0.0, 0.0, 0.0]])-4.0*t2*np.array([[0.0, 0.0, 0.0],[0.0, 0.0, fkk],[0.0, gkk, 0.0]])+dim*np.array([[0.0, -fdim, -gdim],[fdim, 0.0, 0.0],[ gdim, 0.0, 0.0]])
                Ek = -2.0*t1*np.array([[0.0, fk, gk],[fk, 0.0, 0.0],[gk, 0.0, 0.0]])-4.0*t2*np.array([[0.0, 0.0, 0.0],[0.0, 0.0, fkk],[0.0, gkk, 0.0]])+dim*np.array([[0.0, -fdim, 0.0],[fdim, 0.0, 0.0],[0.0, 0.0, 0.0]])
                #vk=LA.eigvals(Ek)
                vk, wk=np.linalg.eigh(Ek)
                #print (wk)
                gapplus[nkx][nky]=vk[0].real
                gapflat[nkx][nky]=vk[1].real
                gapminus[nkx][nky]=vk[2].real
                for nbasis in range(0, basis):
                    wavefn_0[nkx][nky][nbasis]= np.dot(wk[nbasis,0],wk[nbasis,0]).real
                    wavefn_1[nkx][nky][nbasis]= np.dot(wk[nbasis,1],wk[nbasis,1]).real
                    wavefn_2[nkx][nky][nbasis]= np.dot(wk[nbasis,2],wk[nbasis,2]).real
                    g.write("%5.2f %5.2f %5.2f\n" % (wavefn_0[nkx][nky][nbasis],wavefn_1[nkx][nky][nbasis],wavefn_2[nkx][nky][nbasis]))
                f.write("%5.2f %5.2f %5.2f\n" % (gapplus[nkx][nky],gapflat[nkx][nky],gapminus[nkx][nky]))
                
        f.close()
        np.savetxt('kvecs.dat',np.c_[kx_m,ky_m],fmt='%1.4f %1.4f',header="# kx, ky",comments='')

        f = open('dos.dat','w')
        n=0.0
        n1=0.0
        dosbasis=np.zeros((basis,totomega));
        for nbasis in range(0,basis):
            for Nomega in range(1, totomega):
                omega=10.0/totomega*(Nomega-totomega/2)
                dosplus=0.0
                dosflat =0.0
                dosminus=0.0
                for nkx in range (0,Nk):
                    for nky in range(0,Nk):	
                            
                        z=omega+0.02j
                        A1=-(1.0/pi)*(wavefn_0[nkx][nky][nbasis]/(z-gapplus[nkx][nky])).imag
                        A2=-(1.0/pi)*(wavefn_1[nkx][nky][nbasis]/(z-gapflat[nkx][nky])).imag
                        A3=-(1.0/pi)*(wavefn_2[nkx][nky][nbasis]/(z-gapminus[nkx][nky])).imag
                           
                        dosplus=dosplus+A1
                        dosflat=dosflat+A2
                        dosminus=dosminus+A3

                dosbasis[nbasis,Nomega]=dosflat+dosminus+dosplus
                if(nbasis==0):   
                    f.write("%6.2f %6.2f %6.2f\n" %(nbasis,omega, dosbasis[nbasis,Nomega]/(Nk*Nk)))
                
                
                
             
                
        
        


        entropy = open('entropy.dat','w')
        for i in range(1,60):
                temp=(1.0/10)*i
                sum_entropy=0.0
                for Nomega in range(1,totomega):
                        omega=12.0/totomega*(Nomega-totomega/2)
                        dosplus=0.0
                        dosflat =0.0
                        dosminus=0.0
                        for nkx in range (0,Nk):
                                for nky in range(0,Nk):
                                        z=omega+0.02j
                                        A1=-(1.0/pi)*(1.0/(z-gapplus[nkx][nky])).imag
                                        dosplus=dosplus+A1
                                        A2=-(1.0/pi)*(1.0/(z-gapflat[nkx][nky])).imag
                                        dosflat=dosflat+A2
                                        A3=-(1.0/pi)*(1.0/(z-gapminus[nkx][nky])).imag
                                        dosminus=dosminus+A3
                        s=fermi(omega,0.0,temp)
                        sum_entropy=sum_entropy-2.0*(s*log(s))*(dosplus/(Nk*Nk)+dosflat/(Nk*Nk)+dosminus/(Nk*Nk))
                sum_entropy_tot=sum_entropy*(12.0/totomega)
                entropy.write("%6.2f %6.2f\n" % (temp, sum_entropy_tot))
        entropy.close()


#---------------------------------------------------------------------------------------------
def fermi(energy , Ef, temp):
    fermi=1.0/(exp((energy-Ef)/temp)+1.0)
    return float(fermi.real)
#----------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
def onedim(t1,t2,a,Nk,mass,dim,totomega,square=True):
    print (t1,t2,a,Nk)
    
    gap1=np.zeros(Nk)
    gap2=np.zeros(Nk)
    gap3=np.zeros(Nk)
    gap4=np.zeros(Nk)
    f = open("oned.dat", 'w')
    V=mass

    alpha=0.0
    kx_m=[]

    for nkx in range(0, Nk):
    	kx = (2.0*pi)/(Nk-1)*nkx
    	kx_m.append(kx)
    	fk = 1.0+exp(-1j*kx*a)
        
    	gk = 1.0
     	#Ek = np.array([[0.0,(1+alpha)*t1*fk,-(1+alpha)*t1*gk,0],[(1+alpha)*t1*np.conjugate(fk),0.0,0,(1-alpha)*t1*gk],[-(1+alpha)*t1*np.conjugate(gk),0,0.0,(1-alpha)*t1*fk],[0,(1-alpha)*t1*np.conjugate(gk),(1-alpha)*t1*np.conjugate(fk),0.0]])
    	Ek = np.array([[V,(1+alpha)*t1*fk,(1+alpha)*t1*gk,0],[(1+alpha)*t1*np.conjugate(fk),0.0,0,(1-alpha)*t1*gk],[(1+alpha)*t1*np.conjugate(gk),0,0.0,(1-alpha)*t1*fk],[0,(1-alpha)*t1*np.conjugate(gk),(1-alpha)*t1*np.conjugate(fk),0.0]])
    	vk=LA.eigvals(Ek)
    	vk=sorted(vk)

    	gap1[nkx]=vk[0].real
    	gap2[nkx]=vk[1].real
    	gap3[nkx]=vk[2].real
    	gap4[nkx]=vk[3].real

    	f.write("%5.2f %5.2f %5.2f %5.2f\n" % (gap1[nkx],gap2[nkx],gap3[nkx],gap4[nkx]))
    f.close()
    
    np.savetxt('kvecs.dat',np.c_[kx_m],fmt='%1.4f',header="# kx",comments='')


    dos1=np.zeros(totomega)
    dos2=np.zeros(totomega)
    dos3=np.zeros(totomega)
    dos4=np.zeros(totomega)
    f = open('dos.dat','w')
    for Nomega in range(1,totomega):
    	omega=20.0/totomega*(Nomega-totomega/2)
    	dossum1 =0.0
    	dossum2 =0.0
    	dossum3 =0.0
    	dossum4 =0.0

    	for nkx in range (0,Nk):
    		z=omega+0.04j
    		A1=-(1.0/pi)*(1.0/(z-gap1[nkx])).imag
    		dossum1=dossum1+A1
    		A2=-(1.0/pi)*(1.0/(z-gap2[nkx])).imag
    		dossum2=dossum2+A2
    		A3=-(1.0/pi)*(1.0/(z-gap3[nkx])).imag
    		dossum3=dossum3+A3
    		A4=-(1.0/pi)*(1.0/(z-gap4[nkx])).imag
    		dossum4=dossum4+A4
    	dos1[Nomega]=dossum1
    	dos2[Nomega]=dossum2
    	dos3[Nomega]=dossum3
    	dos4[Nomega]=dossum4

    	f.write("%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n" %(omega,dos1[Nomega]/(Nk),dos2[Nomega]/(Nk),dos3[Nomega]/(Nk),dos4[Nomega]/(Nk), dos1[Nomega]/(Nk)+dos2[Nomega]/(Nk)+dos3[Nomega]/(Nk)+dos4[Nomega]/(Nk)))


    h = open("one_gamma_x_m_gamma.dat", 'w')

    for nkx in range(0, Nk):
    	kx = (2*pi)/(Nk)*nkx
    	kx_m.append(kx)
    	h.write("%5.2f %5.2f %5.2f %5.2f %5.2f\n" % (kx,gap1[nkx],gap2[nkx],gap3[nkx],gap4[nkx]))
      			
    h.close()




#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
