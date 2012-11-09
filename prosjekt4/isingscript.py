import os, glob, numpy as np
import matplotlib.pyplot as mpl

N = [10]
j=0
counter = [0]*len(N)
for i in N:
	print "-------HER!----------"
	os.system('./kalle %d'%i)
	outfile = open('collected_results_ising_n_%d.txt'%i,'w')
	counter[j] = i
	kake = []
	for dings in sorted(glob.glob('isingresults_n%d*.txt'%i)):
		noe = open(dings,'r')
		kake.append(noe.read())
		outfile.write(kake[-1])
		os.remove(dings)
	outfile.close()
	j+=1

mpl.figure(1)
n=0
for somefile in sorted(glob.glob('collected_results_ising*.txt')):
	print somefile
	infile = np.loadtxt(somefile)
	mpl.plot(infile[:,-1],infile[:,0],label='%d by %d' %(counter[n],counter[n]))
	mpl.hold('on')
	mpl.xlabel('temperature in units of kT/J')
	mpl.ylabel('average energy per particle')
	#mpl.figlegend((line1),'%d by %d' %(counter[n],counter[n]),'upper left')
	#n +=1
mpl.legend(loc=2)

#mpl.savefig('filename.extension')

mpl.figure(2)
n=0
for somefile in sorted(glob.glob('collected_results_ising*.txt')):
	infile = np.loadtxt(somefile)
	mpl.plot(infile[:,-1],infile[:,1],label='%d by %d' %(counter[n],counter[n]))
	mpl.hold('on')
	mpl.xlabel('temperature in units of kT/J')
	mpl.ylabel('heat capacity per particle')
	
mpl.legend(loc=2)
#mpl.savefig('filename.extension')
mpl.figure(3)
n=0

for somefile in sorted(glob.glob('collected_results_ising*.txt')):
	infile = np.loadtxt(somefile)
	mpl.plot(infile[:,-1],infile[:,2],label='%d by %d' %(counter[n],counter[n]))
	mpl.hold('on')
	mpl.xlabel('temperature in units of kT/J')
	mpl.ylabel('average magnetization per particle')

mpl.legend(loc=1)
#mpl.savefig('filename.extension')
mpl.figure(4)
n=0

for somefile in sorted(glob.glob('collected_results_ising*.txt')):
	infile = np.loadtxt(somefile)
	mpl.plot(infile[:,-1],infile[:,1],label='%d by %d' %(counter[n],counter[n]))
	mpl.hold('on')
	mpl.xlabel('temperature in units of kT/J')
	mpl.ylabel('magnetic suceptibility per particle')

mpl.legend(loc=2)
#mpl.savefig('filename.extension')
mpl.show()
