GNUBND
=========

* Function :
    * Plot bands from bnds.xxx file and syml.xxx file (optional)

* Usage  
    - Type 'gnubnd.py xxx' (where xxx is the last name of bnds file)
    - Type 'gnubnd.py xxx --color for band with color weight

- Output : 
    - bnds.gnu file which is a gnuplot script
    - BNDS.DAT, FERMI.DAT and BNDS2.DAT(spin pol only)
    - after running 'gnuplot bnds.gnu' for postscript output
               'G' will be displayed as symbol Gamma

- Requirement : 
    - Python 2.7 or newer (but should work for all python 2.x)
    - Numpy - numerical library for python
    - To install numpy in ubuntu just use apt-get
        - for mac, download free packages from enthought.com
        - for windows (cygwin), numpy will be available to select when install cygwin
    - Gnuplot
