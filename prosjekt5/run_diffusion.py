import os,argparse,glob,numpy as np
import matplotlib.pyplot as mpl
#from mayavi import mlab
#from mayavi.api import OffScreenEngine

parser = argparse.ArgumentParser()
parser.add_argument("-run",action="store_true", help="run the diffusion.cpp excecutable file")
parser.add_argument("-compile",action="store_true", help="complie the diffusion project")
parser.add_argument("-tofile", action="store_true", help="write results to file (for plotting)")
parser.add_argument("-spacing", type=int, action="store",dest="spacing",default=10, help="write every spacing'th result to file")
parser.add_argument("-removefiles",action="store_true", help="remove the resultfiles")
parser.add_argument("-FE1D", action="store_true",help="Run the Forward Euler discretization in time in 1D")
parser.add_argument("-BE1D", action="store_true",help="Run the Backward Euler discretization in time in 1D")
parser.add_argument("-CN1D", action="store_true",help="Run the Crank Nicolson discretization in time in 1D")
parser.add_argument("-FE2D", action="store_true",help="Run the Forward Euler discretization in time in 2D")
parser.add_argument("-LF2D", action="store_true",help="Run the Leap Frog discretization in time in 2D")
parser.add_argument("-nx", type=int, action="store",dest="nx",default=10, help="number of spacial gridpoints is n^2 in 2D")
parser.add_argument("-nt", type=int, action="store",dest="nt",default=10, help="number of points in time (2D)")
args = parser.parse_args()

tofile = 1 if args.tofile else 0

FE1D = 1 if args.FE1D else 0
BE1D = 1 if args.BE1D else 0
CN1D = 1 if args.CN1D else 0
FE2D = 1 if args.FE2D else 0
LF2D = 1 if args.LF2D else 0

if args.compile:
	os.system('g++ -o willy -O3 diffusion.cpp -larmadillo')

if args.run:
		os.system('./willy %d %d %d %d %d %d %d %d %d' %(tofile,args.spacing,FE1D,BE1D,CN1D,FE2D,\
															LF2D,args.nx,args.nt))



fig = mpl.figure(figsize=(5,5))
ax = fig.add_subplot(111)
for files in sorted(glob.glob('results_BE*.txt')):
	u = np.loadtxt(files)
	#ax.cla()
	#ax.imshow(u)
	mpl.plot(u)
	#picname = files.split('.')
	#picname[1] += '.png'
	#fig.savefig(picname[1])
	mpl.show()
	if args.removefiles:
		os.remove(files)
'''
os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 \
	-ovc lavc -lavcopts vcodec=wmv2 -oac copy -o animation.mpg")
'''