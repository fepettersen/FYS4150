import os,argparse,glob,numpy as np
import matplotlib.pyplot as mpl


parser = argparse.ArgumentParser()
parser.add_argument("-run",action="store_true", help="run the diffusion.cpp excecutable file")
parser.add_argument("-compile",action="store_true", help="complie the diffusion project")
parser.add_argument("-tofile", action="store_true", help="write results to file (for plotting)")
parser.add_argument("-spacing", type=int, action="store",dest="spacing",default=10, help="write every spacing'th result to file")
parser.add_argument("-removefiles",action="store_true", help="remove the resultfiles")
parser.add_argument("-FE", action="store_true",help="Run the Forward Euler discretization in time in 2D")
parser.add_argument("-LF", action="store_true",help="Run the Leap Frog discretization in time in 2D")
parser.add_argument("-nx", type=int, action="store",dest="nx",default=10, help="number of spacial gridpoints is n^2 in 2D")
parser.add_argument("-nt", type=int, action="store",dest="nt",default=10, help="number of points in time (2D)")
args = parser.parse_args()

tofile = 1 if args.tofile else 0
FE = 1 if args.FE else 0
LF = 1 if args.LF else 0

if args.compile:
	os.system('g++ -o willy -O3 diffusion.cpp -larmadillo')

if args.run:
		os.system('./willy %d %d %d %d %d %d' %(tofile,args.spacing,args.FE,args.LF,args.nx,args.nt))

for files in sorted(glob.glob('results_FE*.txt')):
	u = np.loadtxt(files)
	mpl.plot(u)
	mpl.show()
	if args.removefiles:
		os.remove(files)
