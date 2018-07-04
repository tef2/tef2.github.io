#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""
 This file is part of TEF.

    TEF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TEF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TEF.  If not, see <http://www.gnu.org/licenses/>.


(C) 2018 Stratmann D, Pathmanathan JS, Postic G, Rey J, Chomilier J
Contact: dirk.stratmann@sorbonne-universite.fr
"""


#############################################################   LIBRAIRIES   ###########################################################################

import os,sys,string,shutil,math,time
from Bio.PDB import *
from Bio.PDB.NACCESS import *
from decimal import *
import copy
import tef_cleanPDB
import tef_assign1
import tef_show

BUILD_NO = "21"

#SCRIPT_PATH = os.path.dirname(os.path.abspath(argv[0]))
TEF_BIN_PATH = "./src"
NACCESS_BIN_PATH = "/usr/local/naccess"
#NACCESS_LIB_PATH = "/usr/local/naccess"
PYMOL_BIN = "/usr/bin/pymol"
JMOL_PATH = ".."
#JMOL_STYLE = "cartoon"
JMOL_STYLE = "rockets"
#NACCESS_PATH = "/nfs/freetown/user/stratmann/RPBS/naccess" 
#TEF_SEARCH_PATH = "/nfs/freetown/user/stratmann/RPBS/tef_search"


JAVA_SCRIPT_SHOW_HIDE = "\
<script type=\"text/javascript\"><!--\n\
/* Script by: www.jtricks.com\n\
 * Version: 20090221\n\
 * Latest version:\n\
 * www.jtricks.com/javascript/blocks/showinghiding.html\n\
 */\n\
function showPageElement(what)\n\
{\n\
    var obj = typeof what == 'object'\n\
        ? what : document.getElementById(what);\n\
\n\
    obj.style.display = 'block';\n\
    return false;\n\
}\n\
\n\
function hidePageElement(what)\n\
{\n\
    var obj = typeof what == 'object'\n\
        ? what : document.getElementById(what);\n\
\n\
    obj.style.display = 'none';\n\
    return false;\n\
}\n\
\n\
function togglePageElementVisibility(what)\n\
{\n\
    var obj = typeof what == 'object'\n\
        ? what : document.getElementById(what);\n\
\n\
    if (obj.style.display == 'none')\n\
        obj.style.display = 'block';\n\
    else\n\
        obj.style.display = 'none';\n\
    return false;\n\
}\n\
\n\
//--></script>\n\n"


JAVA_SCRIPT_CHECKALL = "\
<script type=\"text/javascript\"><!--\n\
// by Nannette Thacker\n\
// http://www.shiningstar.net\n\
// This script checks and unchecks boxes on a form\n\
// Checks and unchecks unlimited number in the group...\n\
// Pass the Checkbox group name...\n\
// call buttons as so:\n\
// <input type=button name=\"CheckAll\"   value=\"Check All\"\n\
	//onClick=\"checkAll(document.myform.list)\">\n\
// <input type=button name=\"UnCheckAll\" value=\"Uncheck All\"\n\
	//onClick=\"uncheckAll(document.myform.list)\">\n\
// -->\n\
\n\
<!-- Begin\n\
function checkAll(field)\n\
{\n\
for (i = 0; i < field.length; i++)\n\
	field[i].checked = true ;\n\
}\n\
\n\
function uncheckAll(field)\n\
{\n\
for (i = 0; i < field.length; i++)\n\
	field[i].checked = false ;\n\
}\n\
//  End -->\n\
</script>\n\n"

############################################################## Dictionnaires ###########################################################################

# dictionnaire contennant les scores de tous les solutions
dico_score={}
#dico_score[chain] = score of current solution

# dictionnaire contennant les carbones alpha
Dico_CA={}
#Dico_CA[chain][CA1][0]
#Dico_CA[chain_name][num]=[seq_numb,float(x_coord),float(y_coord),float(z_coord)]

# dictionnaire contennant les tef avec leur debut et fin 
Dico_tef={}
# Dico_tef[chain][num][0][0] = beg
# Dico_tef[chain][num][0][1] = end

Dico_tef_cad={}
# Dico_tef[chain][num][0] = CA-CA distance


# dictionnaire contennant la sequence de chaque chaine de la proteine
Dico_seq={}
# Dico_seq[chain_name].append((Dico_AA[resi_aa],seq_numb))
# Dico_seq[chain][num] = (AAtype, resID)

# dictionnaire des acides aminees
Dico_AA = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'MSE':'M', 'CGU':'E', 'UNK':'X'}

# dictionnaires permettant de repertorier les tef en fonction de leur debut 
Dico_index={}
# Dico_index[chain][num] = resBeg
index={}
# index[chain][resBeg] = [tef1, tef2, ...]
# dictionnaire contenant les solutions
solutions={}
# solutions[chain] = tef_list of current solution

# dictionnaire pour les tef a ajouter
# Dprim={}
# Dprim[level] = tef_list of the domain at level
# level = 1 ... n


############################################################ Autres fonctions utiles ###################################################################


def fileExists(f):
	try:
		fileObj = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists


def fileIsExecutable(f):
	#if fileExists(f):
	if os.access(f, os.X_OK):
		return True
	return False


def count_iterable(i):
	return sum(1 for e in i)


def writeTEFlist(filename, chain):
	f = open(filename, "w")
	firstResID = int(Dico_seq[chain][0][1])
	lastResID = int(Dico_seq[chain][-1][1])
	f.write("%d %d\n" % (firstResID, lastResID))
	for i in range(1, len(Dico_tef[chain])+1):
		f.write("%d %d %f\n" % (Dico_tef[chain][i][0][0], Dico_tef[chain][i][0][1], Dico_tef_cad[chain][i][0]))
	f.close()


def readSolution(filename, chain):
	f = open(filename)
	lines = f.readlines()
	read_score = lines[0]
	read_solution = []
	for l in lines[1:]:
		read_solution.append(int(l) + 1)
	solutions[chain] = read_solution
	dico_score[chain] = read_score
	f.close()





# permet de choisir une chaine de la proteine
class ChainSelect(Select):
	def __init__(self, chainID):
		self.chainID = chainID
	def accept_chain(self, chain):
		if chain.get_id() == self.chainID:
			return 1
		else:
			return 0


# fontion pour calculer la distance entre deux carbones alpha
def distance(ca1,ca2):
	dist = (((ca1[1] - ca2[1]) ** 2) + ((ca1[2] - ca2[2]) ** 2) + ((ca1[3] - ca2[3]) ** 2))
	return dist


def writePDBlines(filename, lines):
	if fileExists(filename):
		print "writePDBlines: PDB file " + filename + " already exists - will not overwrite."
		return
	f = open(filename, "w")
	for l in lines:
		f.write(l)
	f.close()
	
def transformMultiPDB(pdb_filename):
	f = open(pdb_filename)
	lines = f.readlines()
	f.close()
	os.rename(pdb_filename, pdb_filename + ".multi")	
	start = -1
	outfile = ""
	numFiles = 0
	for index in range(0,len(lines)):
		l = lines[index]
		if l.split()[0] == "HEADER":
			if start == -1:
				start = index
				outfile = l.split()[-1]
				if len(outfile) == 4:
					outfile += ".pdb"
				else:
					print "ERROR in transformMultiPDB() - len(outfile) != 4"
					sys.exit(-1)
			elif index > start and len(outfile) > 0:
				writePDBlines(outfile, lines[start:index])		
				numFiles += 1
				start = index
			else:
				print "ERROR in transformMultiPDB()"
				sys.exit(-1)
		elif l.split()[0] == "END":
			if start != -1 and index > start:
				writePDBlines(outfile, lines[start:index+1])			
				numFiles += 1
				start = -1
			else:
				print "ERROR in transformMultiPDB() - start == -1 or index > start"
				#sys.exit(-1)
	if numFiles <= 1: # no or only one PDB file written, use original file
		os.rename(pdb_filename + ".multi", pdb_filename)
	return numFiles
	

		

###################################################################### command read ##################################################################################

def main(argv):
	global dico_score, Dico_AA, Dico_CA, Dico_index, dico_score, Dico_seq, Dico_tef, Dico_tef_cad
	#global TEF_BIN_PATH, TEF_SEARCH_PATH, NACCESS_BIN_PATH, NACCESS_LIB_PATH, PYMOL_BIN
	global index, solutions
		
	print ""
	print "TEF2    (build "+ BUILD_NO + ")"
	print "---------------------------------------------------------------------"
	print "This tool decomposes protein 3D structures into sub-domain fragments,"
	print "called 'tightened end fragments' (TEF) or 'closed loops'."
	print "---------------------------------------------------------------------"
	print "(C) 2018 Stratmann D, Pathmanathan JS, Postic G, Rey J, Chomilier J"
	print "Contact: dirk.stratmann@sorbonne-universite.fr"
	print ""
	print ""
	
	
	# verification de la ligne de commande tapee par l'utilisateur
	helpPrint = False
	if len(argv) == 1:
		helpPrint = True
	elif len(argv) > 1:
		if argv[1] == "--help" or argv[1] == "-h" or argv[1] == "-help":
			helpPrint = True
	
	if helpPrint:
		print "Usage: ./TEF2 (-pdbFile FILENAME | -pdb DIRNAME) [-out DIRNAME] [-subdir] [-d CLOSURE]"
                print "              [-min MIN_LENGTH] [-max MAX_LENGTH] [-gw GAP_WEIGHT] [-dw CLOSURE_WEIGHT]"
                print "              [-tw FRAGMENTATION_WEIGHT] [-overlap MAX_OVERLAP] [-gap MAX_GAP] [-pymol]"
                print "              [-debug] [-label RUN_LABEL] [-help]"
                print ""
                print ""
		print "Options:"
		print "    -pdbFile    string    Input PDB file, which must be located in the current directory"
		print "    -pdb        string    Directory with the PDB files to process (default: current directory)"
		print "    -out        string    Create an output directory (which must no already exist)"
		print "    -subdir               Create a sub-directory in the output directory for each PDB file"
		print "    -d          float     Maximum closed loop ends distance (default: 10.0 Å)"
		print "    -min        int       Minimum closed loop length (default: 10 residues)"
		print "    -max        int       Maximum closed loop length (default: 100 residues)"
		#print "          -version 2              old TEF 1.0 or new TEF 2.0 program (1 or 2) (default: 2)"
		#print " >>>>>>>> Only for TEF 2.0 (-version 2) : "
		print "    -gw         float     Gap weight in the score (default: 1.0)"
		print "    -dw         float     Cα-Cα distance weight in the score (default: 1.0)"
		print "    -tw         float     Fragmentation weight in the score (default: 0.0)"
		print "    -overlap    int       Maximum overlap between two closed loops (default: 2 residues)"		
		print "    -gap        int       Maximum gap between two closed loops (default: 100 residues)"
		#print "                       -all            creating one file with all solutions (y or n) (default: y)"
		print "    -pymol                Run PyMOL to generate a PNG file"
		print "    -naccess              Use NACCESS"
		print "    -a          float     Maximum solvent accessibility (default: 25.0%)"
		print "    -label      string    Add run label to the output files"
		print "    -debug                Show debug information"
		print "    -help                 This help"
                print ""
                print ""
                print "Examples:"
                print "    ./TEF2 -pdbFile 1pyo.pdb -out myOUTPUT1"
                print "    ./TEF2 -pdb input_examples/ -out myOUTPUT2"
                print "    ./TEF2 -pdb input_examples/ -out myOUTPUT3 -subdir"
                print ""
		sys.exit()


	#SCRIPT_PATH = os.path.dirname(os.path.abspath(argv[0]))
	#print "SCRIPT_PATH:"
	#print SCRIPT_PATH
	NACCESS_PATH = NACCESS_BIN_PATH + "/naccess" 
	#NACCESS_OPT = " -r " + NACCESS_LIB_PATH + "/vdw.radii -s " + NACCESS_LIB_PATH + "/standard.data"
	#print "NACCESS_PATH:"
	#print NACCESS_PATH
	TEF_SEARCH_PATH = TEF_BIN_PATH + "/tef_search"
	
	##############################################################  PARAMETRES   ###########################################################################
	
	# chevauchement
	MIN=2
	
	# gap lors de la selection du prochain tef
	MAX=100
	
	# pourcentage d'accessibilite maximum d'un residu de la proteine
	MAX_ACCESS = 25
	
	# distance max entre deux carbones alpha
	dist_ca = 10.0 ** 2
	
	# longueur minimum et maximum des tef
	l_min = 10
	l_max = 100
	
	# poids d'un tef pour le calcul du score
	tef_weight=0.0
	
	# poids des gaps pour le calcul du score
	gap_weight=1.0
	
	# poids des distances ca-ca pour le calcul du score
	dist_weight=1.0
	
	# current working directory, can be changed to separate all works
	dirw=os.getcwd()
	pdb_directory=os.getcwd()
	pdb_filename = ""
	
	# save solution of all PDB in the same file by default no "n"
	allinone = "y"
	
	# option for approach 1
	renumber = "F"
	
	# approach default
	approach = "2"
	
	USE_NACCESS = False
	
	USE_SUBDIR = False
	
	USE_PYMOL = False
	
	SHOW_DEBUG = False
	
	runLabel = ""
	
	#################################################################### MAIN ###############################################################################
	
	for i in range(len(argv)):
		if argv[i] == "-pdbFile":
			pdb_filename = argv[i+1]
		elif argv[i] == "-pdb":
			pdb_directory = argv[i+1]
		elif argv[i] == "-d":
			dist_ca = float(argv[i+1])**2
		elif argv[i] == "-min":
			l_min = int(argv[i+1])
		elif argv[i] == "-max":
			l_max = int(argv[i+1])
		elif argv[i] == "-gap":
			MAX = int(argv[i+1])
		elif argv[i] == "-overlap":
			MIN = int(argv[i+1])
		elif argv[i] == "-tw":
			tef_weight = float(argv[i+1])
		elif argv[i] == "-gw":
			gap_weight = float(argv[i+1])
		elif argv[i] == "-dw":
			dist_weight = float(argv[i+1])
		elif argv[i] == "-a":
			MAX_ACCESS= float(argv[i+1])
		elif argv[i] == "-out":
			dirw = argv[i+1]
		#elif argv[i] == "-all":
		#	allinone = argv[i+1]
		elif argv[i] == "-version":
			approach = argv[i+1]
		elif argv[i] == "-naccess":
			USE_NACCESS = True
		elif argv[i] == "-subdir":
			USE_SUBDIR = True
		elif argv[i] == "-pymol":
			USE_PYMOL = True
		elif argv[i] == "-debug":				
			SHOW_DEBUG = True
		elif argv[i] == "-renum" and string.upper(argv[i+1]) == "T":
			renumber = "T"
		elif argv[i] == "-renum" and string.upper(argv[i+1]) == "F":
			renumber = "F"
		elif argv[i] == "-label":
			runLabel = argv[i+1]
		
	
	#Checks
	if USE_NACCESS and not fileIsExecutable(NACCESS_PATH):
		print "ERROR: cannot find naccess executable NACCESS_PATH = "+NACCESS_PATH
		sys.exit(-1)
	if not fileIsExecutable(TEF_SEARCH_PATH):	
		print "ERROR: cannot find tef_search executable TEF_SEARCH_PATH = "+TEF_SEARCH_PATH
		sys.exit(-1)
			
			
	# to check the running time 
	start_time=time.time()
	
	
	
	
	# new work directory
	
	if dirw != os.getcwd():
		if os.path.isdir(dirw) == True :
			print "The directory "+dirw+" exists, please enter a new directory name.\n"
			#print "New directory name :"
			#dirw=raw_input()
			#os.mkdir(dirw)
			exit(-1)
		else:
			os.mkdir(dirw)
	
	pdb_names=[]
	numMultiFiles = -1
	# convert a specific pdb file (also multi PDB) to one or more *.pdb files:
	if len(pdb_filename) > 0:
		if fileExists(pdb_filename):
			# ne copie pas le PDB en int�gralit� (avec tous les mod�les):
			#p = PDB.PDBList(pdb_filename)
			#for a in p:
			#	a.out(outName=a.id.strip()+".pdb")
			
			#alternative:
			#print "pdb_filename=", pdb_filename
			numMultiFiles = transformMultiPDB(pdb_filename)
			#creates one PDB-file per PDB entry in the multi PDB file.
			#files are names according to the PDB code: XXXX.pdb
		else:
			print "ERROR: cannot open file " + pdb_filename
			sys.exit(-1)
	
	# listing of all pdb name
	#print "numMultiFiles=", numMultiFiles
	if numMultiFiles == -1 or numMultiFiles > 1:
		filenames=os.listdir(pdb_directory)
		for name in filenames:
			if len(name.split(".")) == 2:
				if name.split(".")[1] == "pdb":
					if fileExists(pdb_directory + '/' + name):
						pdb_names.append(name)
		pdb_names.sort()
	elif len(pdb_filename) > 0:
		#if len(pdb_filename.split("/")) > 0:
		#	pdb_filename = pdb_filename.split("/")[-1]
		pdb_names.append(pdb_filename)
		
	total_pdb=len(pdb_names)
	
	if total_pdb == 0:
		print "ERROR: no PDB file as input."
		sys.exit(-1)
		
	pdb_teste=0
	
	#print pdb_names
	if approach == "1":
		for name in pdb_names:
			print name
			PDB_id = name.split(".")[0]
			if len(runLabel) > 0:
				PDBlabel = runLabel + "_" + PDB_id
			else:
				PDBlabel = PDB_id
			dmy=string.split(time.asctime())
			hms=dmy[3].split(":")
			if USE_SUBDIR:
				dt=PDBlabel+"_"+dmy[2]+"_"+dmy[1]+"_"+dmy[4]+"_"+hms[0]+hms[1]+hms[2]
				os.mkdir(dirw+"/"+dt) # Create a directory were the results will be saved
				tef_assign1.tef_version1(pdb_directory+"/"+name, dirw+"/"+dt, PDBlabel, math.sqrt(dist_ca), l_min, l_max, renumber)
			else:
				dt=""
				tef_assign1.tef_version1(pdb_directory+"/"+name, dirw, PDBlabel, math.sqrt(dist_ca), l_min, l_max, renumber)
	
		# Check the running time	
		end_time=time.time()
		total_time=end_time-start_time
	
		print "Time: "+str(round(total_time,2))+" s \n"
		print "\n"
	elif approach == "2":	
		# all solution in the same file
		if allinone == "y":
			#os.mkdir(dirw+"/all_results_in_one")
			if len(runLabel) > 0:
				aio=open(dirw+"/"+runLabel+"_all_results.tef","w")
				all_cover=open(dirw+"/"+runLabel+"_all_coverage.tef","w")
				all_ca_dist_mean=open(dirw+"/"+runLabel+"_all_ca_dist_mean.tef","w")
				mean_num_tef=open(dirw+"/"+runLabel+"_all_mean_tef.tef","w")
				f_length=open(dirw+"/"+runLabel+"_all_tef_length.tef","w")
			else:
				aio=open(dirw+"/all_results.tef","w")
				all_cover=open(dirw+"/all_coverage.tef","w")
				all_ca_dist_mean=open(dirw+"/all_ca_dist_mean.tef","w")
				mean_num_tef=open(dirw+"/all_mean_tef.tef","w")
				f_length=open(dirw+"/all_tef_length.tef","w")
		else:
			# file with all TEFs length, used to make a plot with R
			f_length=open(dirw+"/tef_length.tef","w")
			
	
		# create a file with all parameters
		if len(runLabel)>0:
			para_file=open(dirw+"/"+runLabel+"_parameters.tef","w")
		else:
			para_file=open(dirw+"/parameters.tef","w")
		para_file.write("TEF 2.0 (build "+BUILD_NO+") parameters:\n\n")
		para_file.write(' '.join(argv)+"\n\n")
		para_file.write("distance ca-ca (-d) = ")
		para_file.write(str(math.sqrt(dist_ca))+"\n\n")
		para_file.write("length min (-min) = ")
		para_file.write(str(l_min)+"\n\n")
		para_file.write("length max (-max) = ")
		para_file.write(str(l_max)+"\n\n")
		para_file.write("gap weight (-gw) = ")
		para_file.write(str(gap_weight)+"\n\n")
		para_file.write("distance ca-ca weight (-dw) = ")
		para_file.write(str(dist_weight)+"\n\n")
		para_file.write("fragmentation weight (-tw) = ")
		para_file.write(str(tef_weight)+"\n\n")
		para_file.write("max overlap (-overlap) = ")
		para_file.write(str(MIN)+"\n\n")
		para_file.write("max gap (-gap) = ")
		para_file.write(str(MAX)+"\n\n")
		para_file.write("USE_PYMOL (-pymol) = %d\n\n" % USE_PYMOL)
		para_file.write("USE_NACCESS (-naccess) = %d\n\n" % USE_NACCESS)
		para_file.write("max residue accessibility (-a) = ")
		para_file.write(str(MAX_ACCESS)+"\n\n")
		para_file.write("USE_SUBDIR (-subdir) = %d\n\n" % USE_SUBDIR)
		para_file.write("SHOW_DEBUG (-debug) = %d\n\n" % SHOW_DEBUG)		
		para_file.close()
	
		# Start of the program with score approach
		for name in pdb_names:
			renumber="F"
			PDB_id = name.split(".")[0]
			if len(runLabel) > 0:
				PDBlabel = runLabel + "_" + PDB_id
			else:
				PDBlabel = PDB_id
			dmy=string.split(time.asctime())
			hms=dmy[3].split(":")
			if USE_SUBDIR:
				dt=PDBlabel+"_"+dmy[2]+"_"+dmy[1]+"_"+dmy[4]+"_"+hms[0]+hms[1]+hms[2]
				os.mkdir(dirw+"/"+dt) # Create a directory were the results will be saved
				pdbDirw = dirw + "/" + dt
			else:
				dt=""
				pdbDirw = dirw
			model="X"
			#print "Number of PDB to check "+str(total_pdb-pdb_teste)+" .\n" 
			#pdb_teste+=1
			print name
			remove_file = False
			if (fileExists(pdb_directory + '/' + name)==0):
				print "file no longer exists - skip to next file"
				continue
			# cleanPDB is used to check if there is a hole in the PDB
			if (tef_cleanPDB.cleanPDB(pdb_directory+ '/' + name, pdbDirw + '/' + name) == -1):
				print "PDB "+PDB_id+" has got a problem, please check the warning message(s) !\n"
				if USE_SUBDIR:
					os.removedirs(dirw+"/"+dt)
				remove_file = True
				continue
			# Dictionaries
			chains=[]
			dico_score={}
			Dico_CA={}
			Dico_tef={}
			Dico_tef_cad={}
			Dico_seq={}
			Dico_index={}
			index={}
			solutions={}
			# Reading of the PDB, if it's a NMR model we take just the first model
			pdb_file=open(pdbDirw+"/"+name,"r")
			pdb_lines=pdb_file.readlines()
			pdb_file.close()
			for i in range(len(pdb_lines)):
				PDB_renum = pdb_lines[i][0:7]
				if PDB_renum == "REMARK":
					renumber="T"
				check_name=pdb_lines[i][0:5]
				if check_name == "MODEL":
					print "MODEL NMR"
					model = "NMR 1"
					break
			if model != "X":
				tmp_pdb=PDBlabel+"_tmp_NMR.pdb"
				tfo=open(pdbDirw+"/"+tmp_pdb,"w")
				for i in range(len(pdb_lines)):
					check_name=pdb_lines[i][0:6]
					if check_name == "ENDMDL":
						break
					else:
						tfo.write(pdb_lines[i])
				tfo.close()
				tfo=open(pdbDirw+"/"+tmp_pdb,"r")
				pdb_lines=tfo.readlines()
				tfo.close()
			# Listing of Chains
			for i in range(len(pdb_lines)):
				Record_name = pdb_lines[i][0:4]
				if Record_name == "ATOM":
					Atome_name = string.split(pdb_lines[i][12:16])[0]
					chain_name=pdb_lines[i][21]
					if chain_name == " ":
						chain_name = "0"
					if Atome_name == "CA" and not chain_name in chains:
						chains.append(chain_name)
			# Dictionaries	
			for chain in chains:
				print "chain=%s"%(chain)
				Dico_CA[chain]=[]
				Dico_tef[chain]={}
				Dico_tef_cad[chain]={}
				Dico_seq[chain]=[]
				Dico_index[chain]=[]
				index[chain]={}
				solutions[chain]=[]
				dico_score[chain]=[]
			# Listing of CA-atoms
			num=0
			for line in range(len(pdb_lines)):
				Record_name = pdb_lines[line][0:4]
				if Record_name == "ATOM":
					Atome_name = string.split(pdb_lines[line][12:16])[0]
					if Atome_name == "CA":
						num=num+1
						resi_aa=pdb_lines[line][17:20]
						chain_name=pdb_lines[line][21]
						if chain_name == " ":
							chain_name = "0"
						seq_numb=string.split(pdb_lines[line][22:26])[0]
						if not (Dico_AA[resi_aa],seq_numb) in Dico_seq[chain_name]:
							Dico_seq[chain_name].append((Dico_AA[resi_aa],seq_numb))
						x_coord=string.split(pdb_lines[line][30:38])[0]
						y_coord=string.split(pdb_lines[line][38:46])[0]
						z_coord=string.split(pdb_lines[line][46:54])[0]
						Dico_CA[chain_name].append([seq_numb,float(x_coord),float(y_coord),float(z_coord)])
			if model != "X":
				pname=tmp_pdb
			else:
				pname=name
			p = PDBParser(PERMISSIVE=1)	
			s = p.get_structure("test",pdbDirw+"/"+pname)
			# Listing of TEF by chain
			for chain in chains:
				l_sum=0
				l_n=0
				# Creating file with all possible TEFs by chain
				f_position=open(pdbDirw+"/"+PDBlabel+"_positions_chain_"+chain+".tef","w")
				f_position.write("start  end  length  distance ca-ca\n")
				f_position.write("\n")
				num=0
				tmpPDB_id = pname.split(".")[0]
				tmpPDBfile = tmpPDB_id+"_tmp_"+chain+".pdb"
				io=PDBIO()
				io.set_structure(s)
				if chain == "0":
					chain = " "
				io.save(tmpPDBfile, ChainSelect(chain))
				if chain == " ":
					chain = "0"
				tmpModel = p.get_structure(tmpPDB_id+"_"+chain, tmpPDBfile)[0]
				tmpChains = tmpModel.get_list()
				if len(tmpChains) != 1:
					print "ERROR: len(tmpChains) != 1"
					exit(-1)
				atomNum = count_iterable(tmpChains[0].get_atoms())
				resNum = len(tmpChains[0])
				if atomNum < resNum*4:
					atomsMissing = True
					print "Atoms are missing - skip solvant accessibility filter"
				else:
					atomsMissing = False
				surfaceResChain = []
				if not atomsMissing and USE_NACCESS:
					#print NACCESS_PATH+NACCESS_OPT
					#naccess_binary must be just the binary path without any options and whitespaces, 
					#otherwise this will raise an error with Popen from subprocess used by NACCESS class
					#naccessObj = NACCESS(tmpModel, naccess_binary = NACCESS_PATH+NACCESS_OPT)
					naccessObj = NACCESS(tmpModel, naccess_binary = NACCESS_PATH)
					for (resID, values) in naccessObj.__iter__():
						if values['main_chain_rel'] <= MAX_ACCESS and values['side_chain_rel'] <= MAX_ACCESS:
							surfaceResChain.append(resID.get_id()[1])
				for posOne in range(0, len(Dico_CA[chain])):
					for posTwo in range(posOne, len(Dico_CA[chain])):
						l_t=int(Dico_CA[chain][posTwo][0]) - int(Dico_CA[chain][posOne][0])+1
						if (l_t > l_max):
							break
						elif (l_t < l_min):
							continue
						elif not atomsMissing and USE_NACCESS:
							if not (int(Dico_CA[chain][posOne][0]) in surfaceResChain) or not (int(Dico_CA[chain][posTwo][0]) in surfaceResChain):
								continue
						d=distance(Dico_CA[chain][posOne],Dico_CA[chain][posTwo])
						if d > dist_ca:
							continue
						d = math.sqrt(d)
						num=num+1
						Dico_tef[chain][num]=[(int(Dico_CA[chain][posOne][0]),int(Dico_CA[chain][posTwo][0]))]
						Dico_tef_cad[chain][num]=[d]
						dummy = len(Dico_CA[chain][posOne][0])
						if dummy ==1:
							space_1="      "
						elif dummy ==2:
							space_1="     "
						elif dummy ==3:
							space_1="    "
						elif dummy==4:
							space_1="   "
						dummy = len(Dico_CA[chain][posTwo][0])
						if dummy ==1:
							space_2="    "
						elif dummy ==2:
							space_2="   "
						elif dummy ==3:
							space_2="  "
						elif dummy ==4:
							space_2=" "
	
						if len(str(l_t))==1:
							space_3="       "
						elif len(str(l_t))==2:
							space_3="      "
						elif len(str(l_t))==3:
							space_3="     "
						elif len(str(l_t))==4:
							space_3="    "
						ddeci=Decimal(str(d))
						d_s=ddeci.quantize(Decimal('.01'), rounding=ROUND_HALF_UP)
						d_p=str(d_s)
						f_position.write(Dico_CA[chain][posOne][0]+space_1+Dico_CA[chain][posTwo][0]+space_2+str(l_t)+space_3+d_p+" Å\n")
						l_sum=l_sum+l_t
						l_n=l_n+1
				#print "CA-CA finished"
				f_position.write("\n")
				if l_n != 0:
					f_position.write("Mean length: "+str(l_sum/l_n)+"\n")
				elif l_n == 0:
					f_position.write("NO TEFs\n")
				f_position.close()
				os.remove(tmpPDBfile)
				# Use of tef_search a C++ script to do the complete search of all possible best combinations
				if len(Dico_tef[chain]) > 0:
					solution_found={}
					solution_found_score={}
					inputFileName = dirw+"/"+"tef_list_" + tmpPDB_id + "_" + chain + ".txt"
					outputFileName = dirw+"/"+"tef_solution_" + tmpPDB_id + "_" + chain + ".txt"
					#print "TEF_search:"
					writeTEFlist(inputFileName, chain)
					tefSearchCommand = TEF_SEARCH_PATH + " %s %s %i %i %f %f %f %f" % (inputFileName, outputFileName, MIN, MAX, math.sqrt(dist_ca), tef_weight, gap_weight, dist_weight)
					os.system(tefSearchCommand)
					readSolution(outputFileName, chain)
					#print "TEF_search finished"
					if not SHOW_DEBUG:
						os.remove(inputFileName)
						os.remove(outputFileName)
					else:
						print tefSearchCommand
			for chain in chains:
				for tef in solutions[chain]:
					#print str(tef)+"\n"
					length=int(Dico_tef[chain][tef][0][1])-int(Dico_tef[chain][tef][0][0])+1
					f_length.write(str(length)+"\n")
			# Creating a pymol and a Jmol script which allow to see the decomposition of the protein by TEFs
			pymol_file=open(pdbDirw+"/"+PDBlabel+".pymol","w")
			jmol_file=open(pdbDirw+"/"+PDBlabel+".jmol","w")
			pymol_file.write("cmd.load('"+name+"')\n")
			jmol_file.write("load "+name+"\n")
			pymol_file.write("cmd.bg_color('white')\n")
			jmol_file.write("color background white\n")
			pymol_file.write("cmd.hide('everything')\n")		
			pymol_file.write("cmd.show('cartoon')\n")
			jmol_file.write("cartoons Only\n")		
			pymol_file.write("cmd.color('gray')\n")
			pymol_file.write("cmd.set('fog', 'off')\n")
			jmol_html_file=open(pdbDirw+"/"+PDBlabel+".html","w")
			jmol_html_file.write("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n")
			jmol_html_file.write("<html xmlns=\"http://www.w3.org/1999/xhtml\" xml:lang=\"fr\" >\n")
			jmol_html_file.write("	<head>")
			jmol_html_file.write("		<title> TEF 2.0 results on "+ PDBlabel +" </title>\n")
			jmol_html_file.write("		<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n")
			jmol_html_file.write("	</head>\n")
			jmol_html_file.write(JAVA_SCRIPT_SHOW_HIDE)
			jmol_html_file.write(JAVA_SCRIPT_CHECKALL)
			jmol_html_file.write("	<body style=\"font-family:monospace;font-size:11px;line-height:18px;\">\n")
			jmol_html_file.write("<h1> TEF 2.0 results on "+ PDBlabel +" </h1>\n")
			jmol_html_file.write("	<SCRIPT type=\"text/javascript\" src=\""+ JMOL_PATH +"/Jmol.js\"></SCRIPT>\n")#path for jmol folder
			jmol_html_file.write("	<table style=\"font-family:monospace;font-size:12px;line-height:18px;\">\n")
			jmol_html_file.write("		<tr>\n")
			jmol_html_file.write("			<td>\n")
			jmol_html_file.write("				<SCRIPT type=\"text/javascript\" language=\"JavaScript\">\n")
			jmol_html_file.write("					jmolInitialize(\""+ JMOL_PATH+"/\");\n")
			jmol_html_file.write("					jmolApplet(500,\"set echo middle center;font echo 18 sanserif bold;color background white;color echo white;echo Loading molecule||Please wait...; delay 0.1; load "+ name +"; wireframe off; spacefill off; "+JMOL_STYLE+" Only; color beige;\");\n")# need the path of PDB
			jmol_html_file.write("				</SCRIPT>\n")
			jmol_html_file.write("			</td>\n")
			jmol_html_file.write("			<td>\n")
			#jmol_html_file.write("				<table style=\"font-family:monospace;font-size:12px;line-height:18px;\">\n")
			#jmol_html_file.write("					<tr>\n")
			#jmol_html_file.write("						<td>\n")			
			#Creating file with the solution for all chains of the protein
			r=open(pdbDirw+"/"+PDBlabel+"_solutions.tef","w")
			for chain in chains:
				t_n=0
				if chain in solutions:
					sol=solutions[chain]	
					seq_our_h=[]
					seq_our_b=[]
					l_seq=int(Dico_seq[chain][0][1])-1
					if l_seq < 0:
						mis=str(l_seq)
					elif l_seq > 0:
						mis="+"+str(l_seq)
					elif l_seq == 0:
						mis=str(l_seq)
					r.write("PDB  " + PDBlabel + "  Model "+model+"  Chain "+chain+"  Renumber "+renumber+"  shift "+mis+"  score: "+str(dico_score[chain])+"\n")
					r.write("SEQ  ")
					for i in range(len(Dico_seq[chain])):
						r.write(Dico_seq[chain][i][0])
					r.write("\n")
					for i in range(0,len(sol),2):
						catch=sol[i]
						beg=int(Dico_tef[chain][catch][0][0])
						end=int(Dico_tef[chain][catch][0][1])
						seq_our_h.append((beg,end))
					for i in range(1,len(sol),2):
						catch=sol[i]
						beg=int(Dico_tef[chain][catch][0][0])
						end=int(Dico_tef[chain][catch][0][1])
						seq_our_b.append((beg,end))
					shd=int(Dico_seq[chain][0][1])
					r.write("TEF  ")
					if seq_our_h==[]:
						r.write("-"*len(Dico_seq[chain]))
					else:
						for i in range(len(seq_our_h)):
							beg=seq_our_h[i][0]
							end=seq_our_h[i][1]
							if i==0:
								r.write("-"*(beg-shd))
							else:
								r.write("-"*(beg-shd-1))
							r.write("T"*((end-beg)+1))
							shd=end
						r.write("-"*(int(Dico_seq[chain][len(Dico_seq[chain])-1][1])-end))
					r.write("\n")
					sbd=int(Dico_seq[chain][0][1])
					r.write("TEF  ")
					if seq_our_b==[]:
						r.write("-"*len(Dico_seq[chain]))
					else:
						for i in range(len(seq_our_b)):
							beg=seq_our_b[i][0]
							end=seq_our_b[i][1]
							if i== 0:
								r.write("-"*(beg-sbd))
							else:
								r.write("-"*(beg-sbd-1))
							r.write("T"*((end-beg)+1))
							sbd=end
						r.write("-"*(int(Dico_seq[chain][len(Dico_seq[chain])-1][1])-end))
					r.write("\n")
				r.write("\n")
				r.write("TEFs positions:\n")
				r.write("\n")
				for i in range(len(sol)):
					catch=sol[i]
					beg=int(Dico_tef[chain][catch][0][0])
					end=int(Dico_tef[chain][catch][0][1])
					dcc=Dico_tef_cad[chain][catch][0]
					if len(str(i+1))==1:
						space_1="   "
					elif len(str(i+1))==2:
						space_1="  "
					elif len(str(i+1))==3:
						space_1=" "
				
					if len(str(beg))==1:
						space_2="   "
					elif len(str(beg))==2:
						space_2="  "
					elif len(str(beg))==3:
						space_2=" "
	
					if len(str(end))==1:
						space_3="    "
					elif len(str(end))==2:
						space_3="   "
					elif len(str(end))==3:
						space_3="  "
			
					dccdeci=Decimal(str(dcc))
					dcc_s=dccdeci.quantize(Decimal('.01'), rounding=ROUND_HALF_UP)
					dcc_p=str(dcc_s)
					r.write("TEF "+str(i+1)+space_1+"BEG: "+str(beg)+space_2+"END: "+str(end)+space_3+"ca-ca distance: "+dcc_p+" Å"+"  LENGTH: "+str(end-beg+1)+"\n")
					t_n=t_n+1
				if t_n == 0:
					r.write("NO TEFs\n")
				r.write("\n")
				length_seq_chain=len(Dico_seq[chain])
				total_tef_length=0
				s_dcc=0
				colorList = tef_show.color_obj(len(sol),0,0)
				colorListJmol = tef_show.color_obj(len(sol),0,1)
				colorListJmolHex = tef_show.color_obj(len(sol),0,2)
				beg=-9999
				end=-9999
				if len(chains) > 1:
					jmol_html_file.write("<br>\n")

					jmol_html_file.write("<button onclick=\"togglePageElementVisibility('" + chain + "')\">Show/Hide chain " + chain + " results</button><br>")
					jmol_html_file.write("<div id=\""+chain+"\" style=\"display:none; border:1px solid black\">\n")
				jmol_html_file.write("<form name=\"form"+chain+"\">\n")
				jmol_html_file.write("<table>\n")
				jmol_html_file.write("<th colspan=\"5\" align=\"center\">List of TEFs (use checkboxes to color  TEFs):</th>\n");
				jmol_html_file.write("<tr><td>&nbsp;</td><td>View TEF&nbsp;&nbsp;</td><td>begin&nbsp;&nbsp;</td><td>end&nbsp;&nbsp;&nbsp;&nbsp;</td><td>length&nbsp;&nbsp;</td><td>TEF-ends distance</td></tr>\n");
				if chain != 0:
					allButtonText = "jmolScript('display *;select :" + chain + ";color beige;"+JMOL_STYLE+" Only;"
				else:
					allButtonText = "jmolScript('display *;select all;color beige;"+JMOL_STYLE+" Only;"
				for i in range(len(sol)):
					jmol_html_file.write("<tr>\n");
					catch=sol[i]
					oldEnd = end					
					beg=int(Dico_tef[chain][catch][0][0])
					end=int(Dico_tef[chain][catch][0][1])
					dcc=Dico_tef_cad[chain][catch][0]
					s_dcc=s_dcc+dcc
					tef_length_chain=end-beg+1
					total_tef_length=total_tef_length+tef_length_chain
					if beg - oldEnd <= 0:
						total_tef_length -= oldEnd - beg + 1
					jmol_html_file.write("<td>\n");
					if chain != "0":
						pymol_file.write("cmd.color('"+colorList[i]+"','/"+PDB_id+"//"+chain+"/"+str(beg)+"-"+str(end)+"')\n")
						jmol_file.write("select "+str(beg)+"-"+str(end)+":"+chain+"\n")
						jmol_html_file.write("<SCRIPT type=\"text/javascript\" language=\"JavaScript\">jmolCheckbox(\"display *;select "+str(beg)+"-"+str(end)+":"+ chain +";color "+ colorListJmol[i] +"\", \"display *;select "+str(beg)+"-"+str(end)+":"+ chain +";color beige;\", \""+str(i+1)+"\") </SCRIPT>\n")
						allButtonText += "select "+str(beg)+"-"+str(end)+":"+ chain +";color "+ colorListJmol[i] + ";"
					else:
						pymol_file.write("cmd.color('"+colorList[i]+"','/"+PDB_id+"///"+str(beg)+"-"+str(end)+"')\n")
						jmol_file.write("select "+str(beg)+"-"+str(end)+"\n")
						jmol_html_file.write("<SCRIPT type=\"text/javascript\" language=\"JavaScript\">jmolCheckbox(\"display *;select "+str(beg)+"-"+str(end)+";color "+ colorListJmol[i] +"\", \"display *;select "+str(beg)+"-"+str(end)+":"+ chain +";color beige;\", \""+str(i+1)+"\") </SCRIPT>\n")
						allButtonText += "select "+str(beg)+"-"+str(end)+";color "+ colorListJmol[i] + ";"
					jmol_html_file.write("</td>\n");
					jmol_html_file.write("<td>\n");
					jmol_html_file.write("<div class=\"color-box\" style=\"background-color: #"+colorListJmolHex[i]+";\">&nbsp;</div>\n");
					jmol_html_file.write("</td>\n");
					dccdeci=Decimal(str(dcc))
					dcc_s=dccdeci.quantize(Decimal('.01'), rounding=ROUND_HALF_UP)
					dcc_p=str(dcc_s)
					jmol_html_file.write("<td>"+ str(beg) +"</td><td>"+ str(end) +"</td><td>"+ str(tef_length_chain) + "</td><td>"+ dcc_p + "&Aring;</td>\n");
					jmol_file.write("color "+colorListJmol[i]+"\n")
					jmol_html_file.write("</tr>\n");
				jmol_html_file.write("</table>\n")
				jmol_html_file.write("</form>\n")
				allButtonText += "');"
				jmol_html_file.write("<p>\n");
				jmol_html_file.write("<input type=\"button\" name=\"show all TEFs\" value=\"show all TEFs\" onClick=\"checkAll(document.form"+chain+");"+allButtonText+"\"><br>\n")
				if chain != "0":
					notAllButtonText = "jmolScript('display *; select :"+chain+"; color beige');"
				else:
					notAllButtonText = "jmolScript('display *; select all; color beige');"
				jmol_html_file.write("<input type=\"button\" name=\"hide all TEFs\" value=\"hide all TEFs\" onClick=\"uncheckAll(document.form"+chain+");"+notAllButtonText+"\"><br>\n")
				coverage_by_tef=(float(total_tef_length)*100.0)/float(length_seq_chain)
				cbt=Decimal(str(coverage_by_tef))
				cbt_s=cbt.quantize(Decimal('.01'), rounding=ROUND_HALF_UP)
				cbt_p=str(cbt_s)
				r.write("coverage by TEF : "+cbt_p+" %\n")
				r.write("\n")
				jmol_html_file.write("coverage by TEF : "+cbt_p+" %&nbsp;&nbsp;&nbsp;\n")
				if len(sol) == 0:
					s_dcc_p="0"
				else:
					s_dcc=s_dcc/float(len(sol))
					s_dccdeci=Decimal(str(s_dcc))
					s_dcc_s=s_dccdeci.quantize(Decimal('.01'), rounding=ROUND_HALF_UP)
					s_dcc_p=str(s_dcc_s)
				r.write("mean distance ca-ca : "+s_dcc_p+" Å\n")
				jmol_html_file.write("average TEF-ends distance : "+s_dcc_p+" &Aring;&nbsp;&nbsp;&nbsp;\n")
				r.write("\n")
				r.write("sequence length : "+str(length_seq_chain)+"\n")
				jmol_html_file.write("sequence length : "+str(length_seq_chain)+" a.a.<br>\n")
				jmol_html_file.write("</p>\n")
				if len(chains) > 1:
					jmol_html_file.write("</div>\n")
				r.write("\n")
				r.write("\n")
				# Cearting file with solutions of all proteins 
				if allinone == "y":
					all_cover.write(cbt_p+"\n")
					all_ca_dist_mean.write(s_dcc_p+"\n")
					mean_num_tef.write(str(len(sol))+"\n")
					sol=solutions[chain]	
					seq_our_h=[]
					seq_our_b=[]
					l_seq=int(Dico_seq[chain][0][1])-1
					if l_seq < 0:
						mis=str(l_seq)
					elif l_seq > 0:
						mis="+"+str(l_seq)
					elif l_seq == 0:
						mis=str(l_seq)
					aio.write("PDB  "+name+"  Model "+model+"  Chain "+chain+"  Renumber "+renumber+"  shift "+mis+"  score: "+str(dico_score[chain]))#+"\n")
					if dico_score[chain] == []:
						aio.write("\n")
					aio.write("SEQ  ")
					for i in range(len(Dico_seq[chain])):
						aio.write(Dico_seq[chain][i][0])
					aio.write("\n")
					for i in range(0,len(sol),2):
						catch=sol[i]
						beg=int(Dico_tef[chain][catch][0][0])
						end=int(Dico_tef[chain][catch][0][1])
						seq_our_h.append((beg,end))
					for i in range(1,len(sol),2):
						catch=sol[i]
						beg=int(Dico_tef[chain][catch][0][0])
						end=int(Dico_tef[chain][catch][0][1])
						seq_our_b.append((beg,end))
					shd=int(Dico_seq[chain][0][1])
					aio.write("TEF  ")
					if seq_our_h==[]:
						aio.write("-"*len(Dico_seq[chain]))
					else:
						for i in range(len(seq_our_h)):
							beg=seq_our_h[i][0]
							end=seq_our_h[i][1]
							if i==0:
								aio.write("-"*(beg-shd))
							else:
								aio.write("-"*(beg-shd-1))
							aio.write("T"*((end-beg)+1))
							shd=end
						aio.write("-"*(int(Dico_seq[chain][len(Dico_seq[chain])-1][1])-end))
					aio.write("\n")
					sbd=int(Dico_seq[chain][0][1])
					aio.write("TEF  ")
					if seq_our_b==[]:
						aio.write("-"*len(Dico_seq[chain]))
					else:
						for i in range(len(seq_our_b)):
							beg=seq_our_b[i][0]
							end=seq_our_b[i][1]
							if i== 0:
								aio.write("-"*(beg-sbd))
							else:
								aio.write("-"*(beg-sbd-1))
							aio.write("T"*((end-beg)+1))
							sbd=end
						aio.write("-"*(int(Dico_seq[chain][len(Dico_seq[chain])-1][1])-end))
					aio.write("\n")
			r.close()
			pymol_file.write("cmd.ray(1000, 1000)\n")
			pymol_file.write("cmd.png('" + PDBlabel + ".png')\n")
			pymol_file.close()
			#jmol_html_file.write("					</td>\n")
			#jmol_html_file.write("				</tr>\n")
			#jmol_html_file.write("			</table>\n")
			jmol_html_file.write("		</td>\n")
			jmol_html_file.write("	</tr>\n")
			jmol_html_file.write("</table>\n")
			jmol_html_file.write("</body>\n")
			jmol_html_file.write("</html>\n") 
			jmol_html_file.close()			
			jmol_file.close()
			currentDirectory = os.getcwd()
			if USE_PYMOL:
				if fileIsExecutable(PYMOL_BIN):
					os.chdir(pdbDirw)
					os.system(PYMOL_BIN + " -qc " + PDBlabel + ".pymol")
					os.chdir(currentDirectory)
			if model != "X":
				os.remove(pdbDirw +"/"+pname)
			#os.removedirs(dirw+"/"+dt)
			#os.remove(dirw +"/"+name)
		if allinone == "y":
			aio.close()
			all_ca_dist_mean.close()
			mean_num_tef.close()
			f_length.close()
			all_cover.close()
			"""
			# Creating R script
			#r_command=open(dirw +"/all_results_in_one/script_R","w")
			r_command=open(dirw +"/all_script_R","w")
			r_command.write("dca=read.table('all_ca_dist_mean.tef')\n")
			r_command.write("mean(dca)\n")
			r_command.write("coverage=read.table('all_coverage.tef')\n")
			r_command.write("mean(coverage)\n")
			r_command.write("TEF=read.table('all_mean_tef.tef')\n")
			r_command.write("mean(TEF)\n")
			r_command.write("L=read.table('all_tef_length.tef')\n")
			r_command.write("pdf('tef_distribution')\n")
			r_command.write("hist(L[,1],breaks=2000,main='Distribution of TEFs length',xlab='TEF-Length',ylab='Number of TEFs')\n")
			r_command.write("dev.off()\n")
			r_command.write("q()\n")
			r_command.write("n\n")
			r_command.write("\n")
			r_command.close()
			#os.chdir(dirw+"/all_results_in_one")
			#os.system("R < all_script_R --no-save")
			"""
		# Check the running time	
		end_time=time.time()
		total_time=end_time-start_time
	
		print "Time: "+str(round(total_time,2))+" s \n"
		print "\n"



if __name__ == "__main__":
	main(sys.argv)
