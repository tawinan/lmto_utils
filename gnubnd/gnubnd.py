#!/usr/bin/env python
'''
gnubnd program
Created on Oct 31, 2012
Function : Plot bands from bnds.xxx file and syml.xxx file (optional)
Usage    : Type 'gnubnd.py xxx' (where xxx is the last name of bnds file) 
	: Type 'gnubnd.py xxx --color for band with color weight
Output   : bnds.gnu file which is a gnuplot script
	     : BNDS.DAT, FERMI.DAT and BNDS2.DAT(spin pol only) 
	     : after running 'gnuplot bnds.gnu' for postscript output 
	       'G' will be displayed as symbol Gamma  
Requirement : Python 2.7 or newer (but should work for all python 2.x)
            : Numpy - numerical library for python 
            : To install numpy in ubuntu just use apt-get
	        :  for mac, download free packages from enthought.com 
            :  for windows (cygwin), numpy will be available to select 
		        when install cygwin
            : Gnuplot
@author: Tawinan Cheiwchanchamnangij (tawinan.ch@gmail.com)
'''

import sys
from os import system
import numpy as np

EV = 13.60569

def singleLineToFloat(line) :
	"""Input: single line of string 
	   Output :a, b = number of items read, array of items (float)"""
	data = line.split()
	l = len(data)
	x = np.zeros(l)
	for i in range(0,l) :
		x[i] = float(data[i])
	return l, x

def readInt(line):
	"""Input: a single line of string
	Output: an integer read"""
	return int(line.split()[0])

def blockToFloat(block) :
	"""Input: multiple lines of strings 
        Output : a,b = number of items read, 1d array of items (float)"""
	x = []
	l = 0
	for m in block :
		tmp, tmp2 = singleLineToFloat(m)
		l += tmp
		x.append(tmp2)
	en = np.zeros(l)
	i = 0
	j = 0
	for m in x :
		j += len(m)
		en[i:j] = m
		i = j
	return l, en

def readSymLine(lines, size, nEL) :
	""""Input : 
        lines = lines of all k points in one line connecting two high 
		symmetry points
        size = number of k-point in this lines
        nEL = number of line that contains energy values in each k-point
	    Output : a, b = number of k points read, 2D array of energy bands 
		in format which first index (row) refer to k-point, second 
		index (column) refer to energy"""
	nk = np.zeros((size, 3))
	if nEL == 1 :
		bands = np.zeros((size,len(lines[1].split())))
		for i in range(size) :
			ii=i*2
			print "i=", i, "ii=",ii
			print "lines", lines[ii]
			print "lines", lines[ii+1]
			tmp, nktmp = singleLineToFloat(lines[ii])
			nk[i] = nktmp[0:3]
			tmp, bands[i] = singleLineToFloat(lines[ii+1])
	else :
		tmp, tmp2 = blockToFloat(lines[1:(nEL+1)])	
		bands = np.zeros((size, tmp))
		for i in range(size) :
			ii = i*(nEL+1)
			ff = (i+1)*(nEL+1)
			tmp, nktmp = singleLineToFloat(lines[ii])
			nk[i] = nktmp[0:3]
			tmp, bands[i] = blockToFloat(lines[(ii+1):ff])
	return nk, bands
	
def findLength(x1, x2) :
	"""Input: a, b = array of 3d coordinate"""
	return np.sqrt(np.sum((x2 - x1)**2))

def readAllLines(lines, EFERMI, nl, eV = True, ef = True, cb=False):
	"""Input: lines = lines of information, normally start from the second
			 line in bnds file
		  EFERMI = fermi energy in Ry read from the first line of bnds
		  nl = number of lines of band information per 1 k-point, 
			calculated from NB read from first line of bnds file
		  optional eV = True if the output is chosen to be eV 
		  optional ef = True if the output is chosen to reference to 
			fermi energy
            optional cb = True for color band plotting
		Output: lists of 
		  nLine = number of symmetry lines read
		  ebot = minimum energy
		  etop = maximum energy
		  EB = lists of 2D arrays contain energy bands info 
			(see readSymLine for detail)
		  QB = lists of 2D arrays contain k-point read from file
		  NQ = lists of number of k-points read in each symmetry line
          CB = lists of 2D arrays contain band weight info"""
	NQ1 = readInt(lines[0])		#Read first number of k points
	EB = []			#Will be used to store Energy bands
	CB = []         #Will be used to store band weight
	QB = []			#Will be used to store k points
	#Will be used to store number of k points of each sym line
	NQ = []			
	nLine = 0 		#Number of symmetry lines
	startL = 1
	ebot = 0.0
	etop = 0.0
    ##### Start reading the remaining information #########
	while NQ1 > 0 :
		NQ.append(NQ1)
		print " NQ = ", NQ1
		if cb :
		    NQ1 = NQ1*2
		endL = startL + NQ1*(nl+1)
		tmp, tmp2 = readSymLine(lines[startL:endL], NQ1, nl)
		if cb :
		    tmp = tmp[::2,:]    #Trim k-point info by a half
		    tmp3 = tmp2[1::2,:] #Separate the weight info
		    tmp2 = tmp2[::2,:]  #Separate the energy info
		    tmp4 = tmp3 >= 0.0  #Eliminate the negative weight to zero
		    tmp3 = tmp3*tmp4
		#Subtracted by fermi energy if chosen to.
		if ef : tmp2 = tmp2 - EFERMI
		#Change the unit if chosen to.	
		if eV : tmp2 = tmp2*EV			
		# Check for maximum and minimum of bands
		if np.amin(tmp2[:,0]) < ebot : ebot = np.amin(tmp2[:,0])
		if np.amax(tmp2[:,-1]) > etop : etop = np.amax(tmp2[:,-1])
		#Append k-points information to list QB
		QB.append(tmp)
		#Append energy bands information to list EB
		EB.append(tmp2)
		#Append band weight information to list CB
		if cb :
			CB.append(tmp3)
		nLine += 1
		#The start of the next sym line is skip by 1 line of NQ1 info
		startL = endL + 1
		#Read next number of k points, if == 0 will exit the while loop
		NQ1 = readInt(lines[endL])	
			
	return nLine, ebot, etop, EB, QB, CB, NQ
    
def colorBandGen(x, emin, emax, sigma, band, weight) :
	"Generate band weight with the gaussian broadening"
	yres = len(x)
	ex = np.linspace(emin, emax, num=yres)
	for i, w in enumerate(weight) :
		if w > 1e-5 :
			x += w*np.exp(-(ex-band[i])**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))	

if __name__ == "__main__" :

	if len(sys.argv) == 1 :
		print " Usage : Type 'gnubnd.py xxx' (xxx is from bnds.xxx) in order to generate band plot."
		print "         Type 'gnubnd.py xxx --color' to generate band with color weight" 
		sys.exit()
    # Options for plotting format	
	pltmode = 1 			#is output method
	print ' Enter output format:'
	print ' 1 = Postscript'
	print ' 5 = X-Windows (X11)'
	pltmode = raw_input()
	pltmode = int(pltmode)

    # Check command line arguments

	if sys.argv[1] == '' :
		symlfile = 'SYML'
		bndsfile = 'BNDS'
	else :
		symlfile = 'syml.' + sys.argv[1]
		bndsfile = 'bnds.' + sys.argv[1]

    # Check if the bnds file has weight on it by checking --color as a 3rd argv
	bcolor = False
	if len(sys.argv) == 3 :
	    if sys.argv[2] == '--color' :
		bcolor = True
        

# Open input files
	findSyml = True
	try :
		fsyml = open(symlfile, 'r')
		#The information of syml file will be stored in syml
		syml = fsyml.readlines()  	
		fsyml.close()
	#The program still work without syml file
	except IOError :			
		findSyml = False

	fbnds = open(bndsfile, 'r')
	#The information of bnds file will be stored in bnds
	bnds = fbnds.readlines()
	fbnds.close()
	
# set plot title

	ptitle = raw_input("enter title:\n")

# Energy in Ry or eV

	ev = raw_input(
		" energies in Rydberg (f) or eV (t) ? (default is Rydberg)\n")
	if ev == 't' : lev = True
	else : lev = False
	
# set energies relative to EF
	
	ef = raw_input(" energies relative to EF (t)? (default is f)\n")
	if ef == 't' : lef = True
	else : lef = False
	
# set bands connected by lines

#	dot = raw_input(" energies connected by lines (t)? (default is t)\n")
#	if dot == 'f' : ldot = False
#	else : ldot = True
	
	
# read to first line of bnds file 
	temp = bnds[0].split()
	NB, EFERM, ALAT = int(temp[0]), float(temp[1]), float(temp[2])
	
	
	if lev : unit = 'eV'
	else : unit = 'Ry'
	
	print ' NB= %d , EFERM= %10.5f %s, ALAT = %10.5f' %(
		NB, EFERM, unit, ALAT)
	#Number of lines that store energy points, 10 is number of 
	#energy point per line	
	nl = NB/10
	if NB%10 != 0 :
		nl +=1
#Read the remaining information
	nLine, ebot, etop, EB, QB, CB, NQ = readAllLines(bnds[1:], EFERM, nl, eV=lev, ef=lef, cb=bcolor)
#	print "qpts", QB[0][0:], QB[0][1:]
	
#check number of spins, there will be a duplication of k point for the 
#spin polarization calculation
	if findLength(QB[0][0,:], QB[0][1,:]) < 1e-10:	
		nsp = 2
		print " ***Spin polarized calculation is detected. Both spins will be plotted.*** "
	else : nsp = 1
	if lef :
		EFERM = 0.0

	print " EBOT = %10.5f %s ETOP = %10.5f %s EFERM = %10.5f %s " %(
		ebot, unit, etop, unit, EFERM, unit) 
	print " NQ = %d  NLINE= %d " %(np.sum(NQ), nLine)
	
	label = ' '
######## Start reading syml file #####
	lst = []
	if findSyml :
        #Delete all comments, empty lines, and final line
		for i, sl in enumerate(syml):
			if sl[0] == '0' or sl[0] == '#' or sl.split() == [] :
				continue
			else: lst.append(sl)
		syml = lst
		if len(syml[0].split()) == 9 :
			label = ''
			for i, sl in enumerate(syml) :
				if i == len(syml) - 1 : 
					ltmp = sl.split()[-2].replace("'",'')
					ltmp = ltmp.replace('"','')
					label += ltmp
					ltmp = syml[i].split()[-1].replace("'",'')
					ltmp = ltmp.replace('"','')
					label += ltmp
					break
				else : 
					ltmp = sl.split()[-2].replace("'",'')
					ltmp = ltmp.replace('"','')
					label += ltmp
			label = label.upper()
			print " The symmetry labels is found in syml file to be:" 
			print " ", label
						 
		else :
			print " Cannot find the syml labels at the end of syml file!"
			print " Please enter the label of %d characters, ex. GAKAG" %(nLine+1)
			label = raw_input()
	else :
		print "**** File " + symlfile + " can't be opened or not found.****"
		print " Please enter the label of %d characters, ex. GAKAG" %(nLine+1)
		label = raw_input()
######## Finish reading syml file #####

######## Format label #################
	label = label.upper()
	if len(label) > (nLine + 1) :
		label = label[0:(nLine + 1)] 		# Cut the label out if it's too long
	elif len(label) < (nLine + 1) :
		for i in range((nLine + 1) - len(label)) :
			label += " "				# Append the label to the correct length
	lb = []
	for l in label :
		lb.append(l)
		
	label = lb
	
	if pltmode == 1 :
		for i in range(len(label)) :
			if label[i] == 'G' : label[i] = '{/Symbol G}'

######## Finish formatting label
			
	emxn = raw_input("Please enter EMIN, EMAX\n")
	emxn = emxn.replace(',', ' ')
	emax = float(emxn.split()[1])
	emin = float(emxn.split()[0])


####### Check for the lowest and highest bands to be plotted ########
	lmin = 0
	ldiff = 10
	lmax = NB - 1
## Check for lowest band
	while (ldiff > 0) :
		for i in EB :
			diff = emin - np.amax(i[:,lmin])
			ldiff = min(diff,ldiff)
		lmin += 1
## Check for highest band
	ldiff = -10
	while ( ldiff < 0) :
		for i in EB :
			diff = emax - np.amin(i[:,lmax])
			ldiff = max(diff, ldiff)
		lmax -= 1
	lmax += 1				#lmax is undercounted by 1
	print " Bands from %d to %d are plotted.\n" %(lmin, lmax)
    
####### Set up color band #######
	if bcolor :
		sigma_txt = raw_input(" Please enter a line width of the color bands in eV (default is 0.2 eV)\n")
		if sigma_txt == '' :
			sigma = 0.2	#in unit of eV
		else :
			sigma = float(sigma_txt)
		print "Please wait, the program is generating the color bands ... "
		yres = 1001
#		sigma = 0.2         #in unit of eV
		xres = np.sum(NQ)
		imcolor = np.zeros((xres, yres))
		i = 0
		for j in range(len(NQ)) :
		    numk = NQ[j]
		    if nsp == 1 :
			for k in range(numk) :
			    colorBandGen(imcolor[i,:], emin, emax, sigma, EB[j][k,(lmin-1):lmax], CB[j][k,(lmin-1):lmax])
			    i += 1
		    else :
			for k in range(0,numk,2) :
			    colorBandGen(imcolor[i,:], emin, emax, sigma, EB[j][k,(lmin-1):lmax], CB[j][k,(lmin-1):lmax])
			    colorBandGen(imcolor[i,:], emin, emax, sigma, EB[j][k+1,(lmin-1):lmax], CB[j][k+1,(lmin-1):lmax])
			    i += 1
		colmax = np.amax(imcolor)
		print " Maximum weight is %10.6f" %colmax
		norm =raw_input(" Enter the normalization factor (Enter 0 to use the max value above.)\n")
		if norm == '0' :
			imcolor = imcolor*100/colmax
		else :
			nfac = float(norm)
			imcolor = imcolor*100/nfac
		   
	
####### Open output files #######

	fbnds  = open('BNDS.DAT', 'w')
	if nsp == 2 : fbnds2  = open('BNDS2.DAT', 'w')
	ffermi = open('FERMI.DAT', 'w')
	fgnu   = open('bnds.gnu', 'w')
	if bcolor :
		fcb = open('WCOLOR.DAT', 'w')

####### Initialize x axis ########
# x axis will be scale to 1 for the first block
################################## 

	sxk = np.zeros(nLine+1)
	sxk[0] = 0.0
	x = 0.0
	scale = findLength(QB[0][0,:], QB[0][-1,:])
##### Write BNDS.DAT
	for i in range(nLine) :
		if nsp == 1 :
			dx = findLength(QB[i][0,:], QB[i][-1,:])/scale/(NQ[i]-1)
			for j in range(NQ[i]) :	
				fbnds.write("%8.4f\t" %x)
				for k in range(lmin-1,lmax) :
					fbnds.write("%8.4f\t" %EB[i][j,k])
				fbnds.write("\n")
				x +=dx
		else :
			dx = findLength(QB[i][0,:], QB[i][-1,:])/scale/((NQ[i]/2)-1)
			for j in range(NQ[i]/2) :	
				fbnds.write("%8.4f\t" %x)
				fbnds2.write("%8.4f\t" %x)
				for k in range(lmin-1,lmax) :
					fbnds.write("%8.4f\t" %EB[i][j*2,k])
					fbnds2.write("%8.4f\t" %EB[i][j*2 + 1,k])
				fbnds.write("\n")
				fbnds2.write("\n")
				x +=dx
		x -= dx
		sxk[i+1] = x
	fbnds.close()
	if nsp == 2 : fbnds2.close()		

##### Write WCOLOR.DAT ##############
	if bcolor :
		ixres = 0
		x = 0.0
		de = (emax - emin)/(yres-1)
		for i in range(nLine) :
		    if nsp == 1 :
			dx = findLength(QB[i][0,:], QB[i][-1,:])/scale/(NQ[i]-1)
			numk = NQ[i]
		    else :
			dx = findLength(QB[i][0,:], QB[i][-1,:])/scale/((NQ[i]/2)-1)
			numk = NQ[i]/2
		    for j in range(numk) :
			for k in range(yres) :
			    fcb.write("%8.4f\t%8.4f\t%8.4f\n" %(x, emin+k*de, imcolor[ixres,k]))
			ixres += 1
			x += dx
			fcb.write("\n")
		    x -= dx
		fcb.close()
            

##### Write FERMI.DAT
	for i in sxk[1:] :
		ffermi.write("%8.4f\t%8.4f\n" %(i, emin))
		ffermi.write("%8.4f\t%8.4f\n" %(i, emax))
		ffermi.write("\n")
	ffermi.write("%8.4f\t%8.4f\n" %(sxk[0], EFERM))
	ffermi.write("%8.4f\t%8.4f\n" %(sxk[-1], EFERM))
	ffermi.close()
##### Write bnds.gnu
	if pltmode == 1 :
		fgnu.write(" set term postscript enh color font 'Times-Roman,24'\n")
		fgnu.write(" set output 'bnds.ps'\n")
	elif pltmode == 5 :
		fgnu.write(" set term x11\n")
	else :
		fgnu.write(" set term postscript enh color\n")
	if bcolor :
		fgnu.write(" set multiplot\n")
	fgnu.write(" set noxzeroaxis\n")
	fgnu.write(" set tics out\n")
	fgnu.write(" set noxtics\n")
	fgnu.write(" set nokey\n")
	fgnu.write(" set xtics ( ")
	for i in range(len(sxk)-1) :
		fgnu.write("' %s ' %8.4f ," %(label[i], sxk[i]))
	fgnu.write("' %s ' %8.4f )\n" %(label[len(sxk)-1], sxk[-1]))
	if bcolor :
		fgnu.write(" set lmargin screen 0.12\n") 
		fgnu.write(" set rmargin screen 0.85\n") 
		fgnu.write(" set bmargin screen 0.10\n") 
		fgnu.write(" set tmargin screen 0.90\n") 
		fgnu.write(" set pm3d map\n") 
		fgnu.write(" set hidden3d\n") 
		fgnu.write(' set palette defined ( 0 "white", .2 "#70FEBD"  ,40 "blue",70 "red", 100 "white")\n') 
	fgnu.write(" set title '%s'\n" %ptitle)
	fgnu.write(" set yrange [%8.4f : %8.4f] \n" %(emin,emax))
	fgnu.write(" set ylabel 'Energy (%s)'\n" %unit)
	if bcolor :
		fgnu.write(" splot 'WCOLOR.DAT'\n") 
		fgnu.write(" unset xtics\n") 
		fgnu.write(" unset ylabel\n") 
		fgnu.write(" unset ytics\n") 
		fgnu.write(" unset title\n") 

	if nsp == 1 :
		fgnu.write(" set style line 1 lt 1 lw 2 lc rgb 'red'\n")
		fgnu.write(""" plot for[i=2:%d] 'BNDS.DAT' using 1:i with line ls 1, \\
		 'FERMI.DAT' using 1:2 with line lc rgb 'black'\n""" %(lmax-lmin+2))
	else :
		fgnu.write(" set style line 1 lt 1 lw 2 lc rgb 'red'\n")
		fgnu.write(" set style line 2 lt 2 lw 2 lc rgb 'green'\n")
		fgnu.write(" set style line 3 lt 1 lw 1 lc rgb 'black'\n")
		fgnu.write(""" plot for[i=2:%d] 'BNDS.DAT' using 1:i with line ls 1, \\
		for[i=2:%d] 'BNDS2.DAT' using 1:i with line ls 2,\\
		'FERMI.DAT' using 1:2 with line ls 3\n""" %(lmax-lmin+2, lmax-lmin+2))
	if bcolor :
		fgnu.write(" unset multiplot\n")
	if pltmode == 5 :
		fgnu.write(" pause -1 'Press any key to continue'\n")
	fgnu.close()
	system("gnuplot -persist bnds.gnu")
	print " Plot are generated.\n\n"

