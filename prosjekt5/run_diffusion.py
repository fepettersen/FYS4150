import os,argparse,glob,numpy as np
import matplotlib.pyplot as mpl


parser = argparse.ArgumentParser()
parser.add_argument("-run",action="store_true", help="run the diffusion.cpp excecutable file")
parser.add_argument("-compile",action="store_true", help="complie the diffusion project")
parser.add_argument("-tofile", action="store_true", help="write results to file (for plotting)")
parser.add_argument("-spacing", type=int, action="store",dest="spacing",default=10, help="write every spacing'th result to file")
parser.add_argument("-removefiles",action="store_true", help="remove the resultfiles")
args = parser.parse_args()

tofile = 1 if args.tofile else 0
if args.compile:
	os.system('g++ -o willy -O3 diffusion.cpp -larmadillo')

if args.run:
		os.system('./willy %d %d' %(tofile,args.spacing))

for files in sorted(glob.glob('results_FE*.txt')):
	u = np.loadtxt(files)
	mpl.plot(u)
	mpl.show()
	if args.removefiles:
		os.remove(files)

'''
parser.add_argument("-plug",action="store_true", help="Verify that the program can reproduce a square plug exactly in 1D")
parser.add_argument("-gauss",action="store_true", help="Choosing gauss initiat condition")
parser.add_argument("-b", type = float, dest="b", help="Damping coeff")
parser.add_argument("-Lx", type = int, dest="Lx", help="Size of area in x direction")
parser.add_argument("-Ly", type = int, dest="Ly", help="Size of area in y direction")
parser.add_argument("-T", type = int, dest="T", help="Number of timesteps")
parser.add_argument("-Nx", type = int, dest="Nx", help="Number of gridpoints in x direction")
parser.add_argument("-Ny", type = int, dest="Ny", help="Number of gridpoints i y direction")
parser.add_argument("-dt", type = float, dest="dt", help="timestep")
parser.add_argument("-movie", action="store_true",help="")
parser.add_argument("-remove", action="store_true",help="remove the picture and text files produced from plotting")
parser.add_argument("-constant", action="store_true",help="verify that the program can produce a constant solution")
'''