import numpy as np
import matplotlib.pyplot as plt
import pylab as py
#import pyfits
import glob
import struct

class Objet:

	def __init__(self,nom):

	
		
		self.ra=[]
		self.dec=[]
		self.fil=[]
		self.jj=[]
		self.flux=[]
		self.fluxerr=[]
		self.mag=[]
		self.magerr=[]


		self.jjG=[]
		self.jjG_flux=[]
		self.fluxG=[]
		self.fluxerrG=[]
		self.magG=[]
		self.magerrG=[]

		self.jjR=[]
		self.jjR_flux=[]
		self.fluxR=[]
		self.fluxerrR=[]
		self.magR=[]
		self.magerrR=[]

		self.jjI=[]
		self.jjI_flux=[]
		self.fluxI=[]
		self.fluxerrI=[]
		self.magI=[]
		self.magerrI=[]

		self.jjZ=[]
		self.jjZ_flux=[]
		self.fluxZ=[]
		self.fluxerrZ=[]
		self.magZ=[]
		self.magerrZ=[]





		lignes = openFichier(nom)
		self.nomFichierSN = nom
		#print nom
	
		ordreMJD = -1
		ordreMAG = -1
		ordreMAGERR = -1
		ordreFLUXERR = -1
		ordreFLUX = -1
		ordreFLT=-1

						
		for i in range(1,len(lignes)):
			#print lignes[i]
			

			while '  ' in lignes[i]:
				lignes[i]=lignes[i].replace('  ',' ')
			self.fil = lignes[i].split(' ')[10]

			tmp = lignes[i].split(' ')[4].split('"')
			if len(tmp)>1:
				self.ra.append(float(lignes[i].split(' ')[4].split('"')[1]))
				self.dec.append(float(lignes[i].split(' ')[7].split('"')[1]))
			else:
				self.ra.append(float(lignes[i].split(' ')[4].split('"')[0]))
				self.dec.append(float(lignes[i].split(' ')[7].split('"')[0]))

			if self.fil=='sdssr':
				tmp = float(lignes[i].split(' ')[0])
				if tmp in self.jjR or "nan" in str(self.fluxR) :
					pass #print float(lignes[i].split(' ')[3]),"     ",self.fluxR[self.jjR.index(tmp)]

				else:
					tmpR = float(lignes[i].split(' ')[3])

					if(tmpR>-1000 and tmpR<1000.0):

						self.jjR.append(tmp)

						self.fluxR.append(tmpR)

						self.fluxerrR.append(float(lignes[i].split(' ')[6]))
                        
			if self.fil=='sdssi':                
				tmp = float(lignes[i].split(' ')[0])
				if tmp in self.jjI or "nan" in str(self.fluxI) :
					pass #print float(lignes[i].split(' ')[3]),"     ",self.fluxR[self.jjR.index(tmp)]

				else:
					tmpI = float(lignes[i].split(' ')[3])

					if(tmpI>-1000 and tmpI<1000.0):

						self.jjI.append(tmp)

						self.fluxI.append(tmpI)

						self.fluxerrI.append(float(lignes[i].split(' ')[6]))
                        
			if self.fil=='sdssg':
				tmp = float(lignes[i].split(' ')[0])
				if tmp in self.jjG or "nan" in str(self.fluxG) :
					pass #print float(lignes[i].split(' ')[3]),"     ",self.fluxR[self.jjR.index(tmp)]

				else:
					tmpG = float(lignes[i].split(' ')[3])

					if(tmpG>-1000 and tmpG<1000.0):

						self.jjG.append(tmp)

						self.fluxG.append(tmpG)
                        
						self.fluxerrG.append(float(lignes[i].split(' ')[6]))
                        
			if self.fil=='sdssz':
				tmp = float(lignes[i].split(' ')[0])
				if tmp in self.jjZ or "nan" in str(self.fluxZ) :
					pass #print float(lignes[i].split(' ')[3]),"     ",self.fluxR[self.jjR.index(tmp)]

				else:
					tmpZ = float(lignes[i].split(' ')[3])

					if(tmpZ>-1000 and tmpZ<1000.0):

						self.jjZ.append(tmp)

						self.fluxZ.append(tmpZ)
                        
						self.fluxerrZ.append(float(lignes[i].split(' ')[6]))
                        
                      

	def testNbPointsFlux(self, n):
		return len(self.jjR)>n


	def transformListFlux(self, tMin, tMax, nbDecoup, filtre):
		interval = (tMax-tMin)/nbDecoup
		tabFlux = np.zeros((nbDecoup))		
		temps = self.getTempsFlux(filtre)
		flux = self.getFlux(filtre)
		if(len(temps)>0):
			if(min(temps)<tMin or max(temps)>tMax):
				print self.nomFichierSN,"ATTENTION : temps bug ??? =====>", min(temps), max(temps)
		else:
			print self.nomFichierSN, " ::", filtre,  "temps vide ", flux


		for i in range(0, len(temps)):
			t = temps[i]
			f = flux[i]
			caseT = int( (t  - tMin)/(tMax - tMin) * nbDecoup )
			if(caseT<0):
				caseT = 0
			elif(caseT>=nbDecoup):
				caseT = nbDecoup-1
			tabFlux[caseT] = f
		return tabFlux

	def getBorne(self,filtre,ratio):
		tabflux=self.getFlux(filtre)
		tabflux=sorted(tabflux)
		return tabflux[int(len(tabflux)*ratio)]


	def getTempsFlux(self, filtre):
		if(filtre == 'r'):
			return self.jjR
		elif(filtre == 'g'):
			return self.jjG
		elif(filtre == 'i'):
			return self.jjI
		elif(filtre == 'z'):
			return self.jjZ
		else:
			print "filtre bug !!!"
			exit(1)


	def getFlux(self, filtre):
		if(filtre == 'r'):
			return self.fluxR
		elif(filtre == 'g'):
			return self.fluxG
		elif(filtre == 'i'):
			return self.fluxI
		elif(filtre == 'z'):
			return self.fluxZ
		else:
			print "filtre bug !!!"
			exit(1)

	def getAllFlux(self):
		return self.fluxR+self.fluxG+self.fluxI+self.fluxZ

	def getAllTempsFlux(self):
		return self.jjR+self.jjI+self.jjG+self.jjZ
            
 
	def getAllError(self):
		return self.fluxerrR+self.fluxerrG+self.fluxerrI+self.fluxerrZ
            
            
            
def openFichier(nom):
	fd = open(nom,'r')
	lignes_tmp = fd.readlines()
	fd.close()
	return lignes_tmp


def recupererData(lignes):
	tab_supernovae=[]
	for l in range(0,len(lignes)):
		ligne=lignes[l].split('\n')[0]
		tab_supernovae.append(Objet(ligne))

	return tab_supernovae