#!/usr/bin/env python
from numpy import array as a,zeros,sqrt,cos,sin,pi
from scipy.misc import comb as bi

class Kidera:
	def __init__(self,sequence):
		self.seq = sequence
		self.N = len(self.seq)
		self.kidera={}
		self.kidera['A']=a([-1.56,-1.67,-0.97,-0.27,-0.93,-0.78,-0.20,-0.08,0.21,-0.48])
		self.kidera['D']=a([0.58,-0.22,-1.58,0.81,-0.92,0.15,-1.52,0.47,0.76,0.70])
		self.kidera['C']=a([0.12,-0.89,0.45,-1.05,-0.71,2.41,1.52,-0.69,1.13,1.10])
		self.kidera['E']=a([-1.45,0.19,-1.61,1.17,-1.31,0.40,0.04,0.38,-0.35,-0.12])
		self.kidera['F']=a([-0.21,0.98,-0.36,-1.43,0.22,-0.81,0.67,1.10,1.71,-0.44])
		self.kidera['I']=a([-0.73,-0.16,1.79,-0.77,-0.54,0.03,-0.83,0.51,0.66,-1.78])
		self.kidera['V']=a([-0.74,-0.71,2.04,-0.40,0.50,-0.81,-1.07,0.06,-0.46,0.65])
		self.kidera['G']=a([1.46,-1.96,-0.23,-0.16,0.10,-0.11,1.32,2.36,-1.66,0.46])
		self.kidera['H']=a([-0.41,0.52,-0.28,0.28,1.61,1.01,-1.85,0.47,1.13,1.63])
		self.kidera['K']=a([-0.34,0.82,-0.23,1.70,1.54,-1.62,1.15,-0.08,-0.48,0.60])
		self.kidera['L']=a([-1.04,0.00,-0.24,-1.10,-0.55,-2.05,0.96,-0.76,0.45,0.93])
		self.kidera['M']=a([-1.40,0.18,-0.42,-0.73,2.00,1.52,0.26,0.11,-1.27,0.27])
		self.kidera['N']=a([1.14,-0.07,-0.12,0.81,0.18,0.37,-0.09,1.23,1.10,-1.73])
		self.kidera['P']=a([2.06,-0.33,-1.15,-0.75,0.88,-0.45,0.30,-2.30,0.74,-0.28])
		self.kidera['Q']=a([-0.47,0.24,0.07,1.10,1.10,0.59,0.84,-0.71,-0.03,-2.33])
		self.kidera['R']=a([0.22,1.27,1.37,1.87,-1.70,0.46,0.92,-0.39,0.23,0.93])
		self.kidera['S']=a([0.81,-1.08,0.16,0.42,-0.21,-0.43,-1.89,-1.15,-0.97,-0.23])
		self.kidera['T']=a([0.26,-0.70,1.21,0.63,-0.10,0.21,0.24,-1.15,-0.56,0.19])
		self.kidera['W']=a([0.30,2.10,-0.72,-1.57,-1.16,0.57,-0.48,-0.40,-2.30,-0.60])
		self.kidera['Y']=a([1.38,1.48,0.80,-0.56,0.00,-0.68,-0.31,1.03,-0.05,0.53])
		self.aminodict = self.kidera.keys()
		
		self.avgPerm = self._avgPermutation()
		a4 = self._a4()
		print a4-self.avgPerm*self.avgPerm
		self.stddev = sqrt(a4-self.avgPerm*self.avgPerm)
	
	def avgPermutation(self):
		return self.avgPerm
	def standardDev(self):
		return self.stddev
			
	def avgFactors(self):
		output=zeros(10)
		l=1.0*len(self.seq)
		for aa in self.seq:
			if aa in self.aminodict:
				output+=self.kidera[aa]
			else:
				print "ERRORRR! ",self.seq,aa
		return output/l	
	def _ksi1(self):
		N = self.N
		aa = self.aminodict
		AAlen = len(aa)
		seq = self.seq
		suma=zeros(10)
		for fidx in range(10):			
			for x in range(AAlen):
				xc = seq.count(aa[x])
				# X^4 * nx(nx-1)*(nx-2)*(nx-3)
				sk = self.kidera[aa[x]][fidx]
				sk*=sk
				sk*=sk
				
				suma[fidx] += sk*xc
		suma /= N
		return suma
		
		
	def _ksi2(self):
		N = self.N
		aa = self.aminodict
		AAlen = len(aa)
		seq = self.seq
		suma=zeros(10)
		for fidx in range(10):
			for x in range(AAlen):
				xc = seq.count(aa[x])
				suma[fidx]+=self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*xc*(xc-1)
				for y in range(AAlen):
					if x!=y:
						yc = seq.count(aa[y])
						suma[fidx]+=self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[y]][fidx]*xc*yc
		suma /= N*(N-1)
		return suma
		
	def _ksi3(self):
		N = self.N
		aa = self.aminodict
		AAlen = len(aa)
		seq = self.seq
		suma=zeros(10)
		for fidx in range(10):
			for x in range(AAlen):
				xc = seq.count(aa[x])
				suma[fidx]+=self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*xc*(xc-1)
				for y in range(AAlen):
					if x!=y:
						yc = seq.count(aa[y])
						suma[fidx]+=self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[y]][fidx]*self.kidera[aa[y]][fidx]*xc*yc
		suma /= N*(N-1)
		return suma
		
	def _ksi4(self):
		N = self.N
		aa = self.aminodict
		AAlen = len(aa)
		seq = self.seq
		suma = zeros(10)
		
		for fidx in range(10):	
			for x in range(AAlen):
				xc = seq.count(aa[x])
				# X^4 * nx(nx-1)*(nx-2)
				vx=self.kidera[aa[x]][fidx]
				suma[fidx] += vx*vx*vx*vx*xc*(xc-1)*(xc-2) 		
				for y in range(AAlen):
					yc = seq.count(aa[y])
					if x!=y and yc>0:
						# 2.0*X^3*Y*nx(nx-1)ny
						suma[fidx]+= 2.0*vx*vx*vx*self.kidera[aa[y]][fidx]*xc*(xc-1)*yc
						# X^2 *Y^2*nx*ny(ny-1)
						suma[fidx]+= self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]*self.kidera[aa[y]][fidx]*self.kidera[aa[y]][fidx]*xc*yc*(yc-1)
						for z in range(AAlen):
							zc = seq.count(aa[z])
							if x!=z and y!=z and zc>0:
								# 2* X*X*Y*Z*nxnynz
								suma[fidx] += vx*vx*self.kidera[aa[y]][fidx]*self.kidera[aa[z]][fidx]*xc*yc*zc
		suma /= N*(N-1)*(N-2)
		return suma
	
	def _ksi5(self):
		N = self.N
		aa = self.aminodict
		AAlen = len(aa)
		seq = self.seq
		suma = zeros(10)
		
		for fidx in range(10):	
			for x in range(AAlen):
				xc = seq.count(aa[x])
				xx = self.kidera[aa[x]][fidx]*self.kidera[aa[x]][fidx]
				suma[fidx] += 4.0*xx*xx*xc*(xc-1)*(xc-2)*(xc-3)
				for y in range(AAlen):
					yc = seq.count(aa[y])
					if  yc>0 and x!=y:
						suma[fidx]+= 48.0*xx*self.kidera[aa[x]][fidx]*self.kidera[aa[y]][fidx]*xc*(xc-1)*yc*(xc-2)
						# kurwa
						suma[fidx]+= 36.0*xx*self.kidera[aa[y]][fidx]*self.kidera[aa[y]][fidx]*xc*yc*(yc-1)*(xc-1)
						
						for z in range(AAlen):
							zc = seq.count(aa[z])
							if  zc>0 and x!=z and y!=z:
								suma[fidx] += 144.0*xx*self.kidera[aa[y]][fidx]*self.kidera[aa[z]][fidx]*xc*yc*zc*(xc-1)
								for p in range(AAlen):
									pc = seq.count(aa[p])
									if pc>0 and p!=z and x!=p and p!=y:
										suma[fidx]+= 24.0*self.kidera[aa[x]][fidx]*self.kidera[aa[y]][fidx]*self.kidera[aa[z]][fidx]*self.kidera[aa[p]][fidx]*xc*yc*zc*pc
		suma /= N*(N-1)*(N-2)*(N-3)
		return suma

	def _S(self, i):
		if i%2!=0:
			return 0.0
		return self.N*bi(i,i/2,1)/(2.0**i)
	def _Sij(self, i,j):
		return self._S(i)*self._S(j) - self._S(i+j)
	def _S211(self):
		return  self.N*(3-self.N)/4.0
	def _S1111(self):
		return  self.N*(-2.25 + 0.75*self.N)
    
	def _a4(self):
		a = zeros(10)
		a = self._S(4)*self._ksi1() + 4.0*self._Sij(3,1)\
		*self._ksi2()+6.0*self._Sij(2,2)*self._ksi3()+12.0*self._S211()*self._ksi4()+6*self._S1111()*self._ksi5()
		return a
		    
	def _avgPermutation(self):
		avg = zeros(10) # avg[0] contains <(a_k^0)^2>, avg over seq permutations of 0-th property
		N = self.N
		Xl = len(self.aminodict)
		aa = self.aminodict
		
		tri = 1.0/((N-1))
		for fct_idx in range(10):
			for X in range(Xl):
				n_x = self.seq.count(aa[X])
				x2 = self.kidera[aa[X]][fct_idx] * self.kidera[aa[X]][fct_idx]
				avg[fct_idx] -= tri*x2 * n_x*(n_x-1) -  x2 * n_x
				for Y in range(Xl):
					if X!=Y:
						n_y = self.seq.count(aa[Y])
						xy = self.kidera[aa[X]][fct_idx] * self.kidera[aa[Y]][fct_idx]
						avg[fct_idx] -= tri*xy*n_x*n_y
		avg*=0.5
		return avg
	def sineFourier(self,k):
		sine=zeros(10)
		for fct_idx in range(10):
			for i in range(self.N):
				sine[fct_idx] += self.kidera[self.seq[i]][fct_idx] * sin(2*pi*i*k/self.N)	
		return sine
	def cosineFourier(self,k):
		cosine=zeros(10)
		for fct_idx in range(10):
			for i in range(self.N):
				cosine[fct_idx] += self.kidera[self.seq[i]][fct_idx] * cos(2*pi*i*k/self.N)	
		return cosine
		
		
		
				
	def ksi5_wolne(self):
		avg = zeros(10) # avg[0] contains <(a_k^0)^2>, avg over seq permutations of 0-th property
		N = len(self.seq)
		Xl = len(self.aminodict)
		aa = self.aminodict
		
		from itertools import permutations,product
		import itertools
		seq0 = list(permutations(self.seq))
		for fct_idx in range(1):
			for seq in seq0:
				for m in range(N):
					for p in range(N):
						if m!=p:
							for n in range(N):
								if n!=p and n!=m:
									for q in range(N):
										if  q!=n and q!=m and q!=p :
											fm = self.kidera[seq[m]][fct_idx]
											fp = self.kidera[seq[p]][fct_idx]
											fn = self.kidera[seq[n]][fct_idx]
											fq = self.kidera[seq[q]][fct_idx]
											avg[fct_idx] += fm*fq*fp*fn*sin(2*pi*m/N)*sin(2*pi*q/N)*sin(2*pi*p/N)*sin(2.0*pi*n/N)

			avg[fct_idx]/=len(seq0)
		return avg
	def ksi4_wolne(self):
		avg = zeros(10) # avg[0] contains <(a_k^0)^2>, avg over seq permutations of 0-th property
		N = len(self.seq)
		Xl = len(self.aminodict)
		
		from itertools import permutations,product
		seq0 =  list(permutations(self.seq))
		for fct_idx in range(1):	
	
			
			for seq in seq0:
				
				for m in range(N):
					for p in range(N):
						for n in range(N):
							if n!=m:
								if n!=p:
									if p!=m:
										fm = self.kidera[seq[m]][fct_idx]
										fp = self.kidera[seq[p]][fct_idx]
										fn = self.kidera[seq[n]][fct_idx]
										avg[fct_idx] += fm*fm*fp*fn*sin(2*pi*m/N)*sin(2*pi*m/N)*sin(2*pi*p/N)*sin(2.0*pi*n/N)
			avg[fct_idx]/=len(seq0)*1.0
				
				
		return avg

	def ksi3_wolne(self):
		avg = zeros(10) # avg[0] contains <(a_k^0)^2>, avg over seq permutations of 0-th property
		N = len(self.seq)
		Xl = len(self.aminodict)
		aa = self.aminodict
		
		for fct_idx in range(1):	
			# for wszystkie permutacje self.seq !!! BUAAHAHAHAH
			from itertools import permutations,product
			seq0 = list(permutations(self.seq))
			for seq in seq0:
				for m in range(N):
					for n in range(N):
						if n!=m:
							fm = self.kidera[seq[m]][fct_idx]
							fn = self.kidera[seq[n]][fct_idx]
							avg[fct_idx] += fm*fm*fn*fn*sin(2*pi*m/N)*sin(2*pi*m/N)*sin(2*pi*n/N)*sin(2.0*pi*n/N)
					
			avg[fct_idx]/=len(seq0)*1.0
				
		return avg
	def ksi1_wolne(self):
		avg = zeros(10) # avg[0] contains <(a_k^0)^2>, avg over seq permutations of 0-th property
		N = len(self.seq)
		Xl = len(self.aminodict)
		aa = self.aminodict
		
		for fct_idx in range(1):	
			# for wszystkie permutacje self.seq !!! BUAAHAHAHAH
			from itertools import permutations,product
			seq0 = list(permutations(self.seq))
			for seq in seq0:
				for m in range(N):
					fm = self.kidera[seq[m]][fct_idx]
					avg[fct_idx] += fm*fm*fm*fm*sin(2*pi*m/N)*sin(2*pi*m/N)*sin(2*pi*m/N)*sin(2.0*pi*m/N)
					
			avg[fct_idx]/=len(seq0)*1.0
				
		return avg

		
	def ksi2_wolne(self):
		avg = zeros(10) # avg[0] contains <(a_k^0)^2>, avg over seq permutations of 0-th property
		N = len(self.seq)
		Xl = len(self.aminodict)
		aa = self.aminodict
		
		for fct_idx in range(1):	
			# for wszystkie permutacje self.seq !!! BUAAHAHAHAH
			from itertools import permutations,product
			seq0 = list(permutations(self.seq))
			for seq in seq0:
				for m in range(N):
					for n in range(N):
						if n!=m:
							fm = self.kidera[seq[m]][fct_idx]
							fn = self.kidera[seq[n]][fct_idx]
							avg[fct_idx] += fm*fm*fm*fn*sin(2*pi*m/N)*sin(2*pi*m/N)*sin(2*pi*m/N)*sin(2.0*pi*n/N)
					
			avg[fct_idx]/=len(seq0)*1.0
				
		return avg




	

T0498 = Kidera("TTYKLILNLKQAKEEAIKELVDAGTAEKYIKLIANAKTVEGVWTLKDEIKTFTVTE")
T0499 = Kidera("TTYKLILNLKQAKEEAIKEAVDAGTAEKYFKLIANAKTVEGVWTYKDEIKTFTVTE")
for k in range(28):
	avgp = T0499.avgPermutation()
	avgstddev = T0499.standardDev()
	co = T0499.sineFourier(k)
	alfa = (co*co - avgp)/avgstddev
	for fi in range(len(alfa)):
		if alfa[fi]>2:
			print "cos",alfa[fi],"f=",fi,"k=",k

