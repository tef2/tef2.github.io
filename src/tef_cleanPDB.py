#!/usr/bin/env python
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



import os,sys
import warnings
from Bio.PDB import *

MIN_DISTANCE = 1
MAX_DISTANCE = 4.2

def fileExists(f):
	try:
		fileObj = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists


def cleanPDB(fileIn, fileOut):		
	p = PDBParser(PERMISSIVE=1)	
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		s = p.get_structure("dummy", fileIn)
	renumber = False
	for model in s:
		for chain in model:
			oldResNo = -9999
			for residue in chain:
				"""
				if residue.is_disordered():
					print "Disordered atoms - skip file %s" % (fileIn)
					return -1
				"""
				resID = residue.get_id()
				hetero_flag = resID[0]
				if hetero_flag != " ":
					continue
				resNo = resID[1]
				if oldResNo != -9999:
					if (oldRes.child_dict.has_key('CA') and residue.child_dict.has_key('CA')):
						ca1=oldRes['CA']
						ca2=residue['CA']								
						distance = ca1 - ca2
						if distance < MIN_DISTANCE or distance > MAX_DISTANCE:
							print "Chain break between resNo = %d and resNo = %d, distance = %2.2f - skip file %s" % (oldResNo, resNo, distance, fileIn)
							return -1
					if resNo != (oldResNo + 1):
							residue.id = (resID[0], oldResNo + 1, resID[2])
							oldResNo += 1
							renumber = True
					else:
						oldResNo = resNo
				else:
					oldResNo = resNo
				oldRes = residue
	io=PDBIO()
	io.set_structure(s)
	if fileExists(fileOut):
		os.system("cp " + fileOut + " " + fileOut + ".org")
	f = open(fileOut, "w")
	if renumber:
		f.write("REMARK PDB has been renumbered\n")
	io.save(f)
	f.close()
	with warnings.catch_warnings():
		warnings.simplefilter("error")
		try:
			p = PDBParser(PERMISSIVE=0)
			s = p.get_structure("dummy", fileOut)
		except:
			print "Error in PDB:\n", sys.exc_info()[1]
			return -1
	return 0


'''#### main
if len(sys.argv) < 5:
	print "Usage: "+sys.argv[0]+" -in InputDirectory -out OutputDirectory"
	sys.exit()


for i in range(len(sys.argv)):
	if sys.argv[i] == "-in":
		input_directory = sys.argv[i+1]
	elif sys.argv[i] == "-out":
		output_directory = sys.argv[i+1]
		
		
pdb_names=os.listdir(input_directory)

numWritten = 0
for name in pdb_names:
	#print input_directory +'/' + name
	sys.stderr.write(input_directory +'/' + name + '\n')
	if (cleanPDB(input_directory + '/' + name, output_directory +'/' + name) == 0): 
		numWritten += 1

print "%d PDB files have been written to %s" % (numWritten, output_directory)'''
