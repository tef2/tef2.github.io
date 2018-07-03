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

#################################################################################################
# Authors: M. Lonquety, N. Papandreou, M. Lamarine, J. Chomilier                                                                                 #
#                                                                                               #
#              Script qui assigne les TEF a partir d'un fichier PDB                             #
#                ./tef_assign.py -pdb fic_pdb -o rep_out                                        #
#                      options : -seuil [float]  distance max entre les extremites              #
#                                -lmin  [int]    longueur minimum d'un TEF                      #
#                                -lmax  [int]    longueur maximum d'un TEF                      #
#                                -renum [F ou T] si on renumerote ou pas les residus            #
#                                -mod   [int]    numero du modele a afficher                    #
#                                -chain [char]   nom de la chaine a afficher                    #
#                                                                                               #
#                                                                                               #
#################################################################################################




import os, sys, string, re, math

dico_aa = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'MSE':'M', 'CGU':'E', 'UNK':'X'}

# fonction qui verifie si des heteroatomes en plein milieu de la sequence doivent etre pris en compte
def verif_modres(ligne, dico_modres):
	verdict = False
	if ligne[0:6] == "HETATM" and ligne[13:15] == "CA":
		if dico_modres.has_key(ligne[17:20]):
			for res in dico_modres[ligne[17:20]]:
				if res[1] == int(ligne[22:26]):
					verdict = True

	return verdict


# fonction d'extraction des coordonnees dans un fichier pdb
def coord_pdb(lignes_pdb):

	modele="1"
	num_res_prec = "NA"
	dico_model={}
	dico_modres={}
	dico_seqres={}
	dico_contigu={}
	dico_missing_res={}
	cpt_res=1
	for li in range(len(lignes_pdb)):
		if lignes_pdb[li][0:5] == "MODEL":
			modele = lignes_pdb[li][10:14]
			dico_model[modele]={}
			dico_contigu[modele]={}
			dico_missing_res[modele]={}
		# il faut prendre en compte le cas ou il y a des heteroatomes en plein milieu de la structure
		# qui sont des residus speciaux et non pas des ions
		elif lignes_pdb[li][0:6] == "MODRES":
			if dico_modres.has_key(lignes_pdb[li][12:15]):
				dico_modres[lignes_pdb[li][12:15]].append([lignes_pdb[li][16], int(lignes_pdb[li][18:22])])
			else:
				dico_modres[lignes_pdb[li][12:15]] = [[lignes_pdb[li][16], int(lignes_pdb[li][18:22])]]
		# on ne s'interesse qu'aux champs ATOM ou sinon aux HETATM specifiques definis par MODRES
		elif (lignes_pdb[li][0:4] == "ATOM" and lignes_pdb[li][13:15] == "CA") or verif_modres(lignes_pdb[li], dico_modres):
			chaine = lignes_pdb[li][21]
			num_res = lignes_pdb[li][22:26]
			if dico_model.has_key(modele):
				# on prend en compte chaque chaine bien entendu
				if dico_model[modele].has_key(chaine):
					# on verifie que les residus sont contigus
					#print num_res_prec+" "+num_res
					if num_res_prec == "NA" or int(num_res_prec) == int(num_res)-1:
						#print " --- > banjo<br>"
						# on rajoute dans la liste les coordonnees x, y, z, le numero du residu et le residu en lui meme
						dico_model[modele][chaine].append([float(lignes_pdb[li][30:38]), float(lignes_pdb[li][38:46]), float(lignes_pdb[li][46:54]), int(lignes_pdb[li][22:26]), lignes_pdb[li][17:20]])
						num_res_prec = num_res
					# si les residus ne sont pas contigus on va chercher a savoir s'il y a reellement un trou dans la structure
					# en calculant la distance separant les 2 residus et en regardant si elle est a 3.8 A +/- 10%
					else:
						# on verifie qu'une forme du residu n'a pas deja ete definie precedemment avant toute chose
						if lignes_pdb[li][16] != " ":
							#print " --- > on zappe<br>"
							continue
						dico_contigu[modele][chaine] = False
						#print "--- > pas bon<br>"
						res_calc_dist = calcul_dist_ext([dico_model[modele][chaine][-1], [float(lignes_pdb[li][30:38]), float(lignes_pdb[li][38:46]), float(lignes_pdb[li][46:54]), int(lignes_pdb[li][22:26]), lignes_pdb[li][17:20]]], 100.0, 2, 1100, renumber)
						test_contigu = res_calc_dist[0]
						#print "<br>%lf &Aring;<br>" %(test_contigu[0][2])
						if test_contigu != [] and test_contigu[0][2] >= 3.4 and test_contigu[0][2] <= 4.2 or lignes_pdb[li][26] != " ":
							dico_model[modele][chaine].append([float(lignes_pdb[li][30:38]), float(lignes_pdb[li][38:46]), float(lignes_pdb[li][46:54]), int(lignes_pdb[li][22:26]), lignes_pdb[li][17:20]])
							num_res_prec = num_res
						else:
							num_res_prec = num_res
							dico_missing_res[modele][chaine] = True
				else:
					dico_contigu[modele][chaine] = True
					dico_model[modele][chaine] = [[float(lignes_pdb[li][30:38]), float(lignes_pdb[li][38:46]), float(lignes_pdb[li][46:54]), int(lignes_pdb[li][22:26]), lignes_pdb[li][17:20]]]
					dico_missing_res[modele][chaine] = False
					num_res_prec = num_res
			else:
				dico_model[modele] = {}
				dico_contigu[modele] = {}
				dico_missing_res[modele] = {}
				dico_contigu[modele][chaine] = True
				dico_model[modele][chaine] = [[float(lignes_pdb[li][30:38]), float(lignes_pdb[li][38:46]), float(lignes_pdb[li][46:54]), int(lignes_pdb[li][22:26]), lignes_pdb[li][17:20]]]
				dico_missing_res[modele][chaine] = False
				num_res_prec = num_res

		# au cas ou, on re initialise num_res_prec qd on change de modele
		if lignes_pdb[li][0:6] == "ENDMDL":
			num_res_prec = "NA"

		# on recupere le header du pdb pour l'afficher sur la page
		#elif lignes_pdb[li][0:6] == "HEADER":
			#print lignes_pdb[li][7:]

		# on recupere la sequence de seqres en cas de trou pour tout avoir
		elif lignes_pdb[li][0:6] == "SEQRES":
			if dico_seqres.has_key(lignes_pdb[li][11]):
				dico_seqres[lignes_pdb[li][11]] += " "+lignes_pdb[li][19:70]
			else:
				dico_seqres[lignes_pdb[li][11]] = lignes_pdb[li][19:70]
	return dico_model, {}, dico_contigu, dico_missing_res

# fonction de calcul de distances entre les extremites de tous les fragments possibles
# on ne retient que ceux qui sont inferieurs a seuil
def calcul_dist_ext(liste_coord, seuil, long_min, long_max, renumber):
	liste_closed_loops=[]
	dico_interactions = {}
	for i in range(len(liste_coord)):
		for j in range(i+1, len(liste_coord)):
			dist = math.sqrt(math.pow((liste_coord[i][0] - liste_coord[j][0]), 2) + math.pow((liste_coord[i][1] - liste_coord[j][1]), 2) + math.pow((liste_coord[i][2] - liste_coord[j][2]), 2))
			#print dist
			if dist <= seuil and (j-i+1) >= long_min and (j-i+1) <= long_max:

				# on stocke les closed loops et en meme temps on note les residus en interaction
				if renumber:
					liste_closed_loops.append([i+1, j+1, dist, code_aa(dico_aa, liste_coord[i][4]), code_aa(dico_aa, liste_coord[j][4])])
					if dico_interactions.has_key(i+1):
						dico_interactions[i+1][0] += 1
					else:
						dico_interactions[i+1] = [1, code_aa(dico_aa, liste_coord[i][4]), i+1]
					if dico_interactions.has_key(j+1):
						dico_interactions[j+1][0] += 1
					else:
						dico_interactions[j+1] = [1, code_aa(dico_aa, liste_coord[j][4]), j+1]
				else:
					liste_closed_loops.append([liste_coord[i][3], liste_coord[j][3], dist, code_aa(dico_aa, liste_coord[i][4]), code_aa(dico_aa, liste_coord[j][4])])
					if dico_interactions.has_key(liste_coord[i][3]):
						dico_interactions[liste_coord[i][3]][0] += 1
					else:
						dico_interactions[liste_coord[i][3]] = [1, code_aa(dico_aa, liste_coord[i][4]), liste_coord[i][3]]
					if dico_interactions.has_key(liste_coord[j][3]):
						dico_interactions[liste_coord[j][3]][0] += 1
					else:
						dico_interactions[liste_coord[j][3]] = [1, code_aa(dico_aa, liste_coord[j][4]), liste_coord[j][3]]
	return liste_closed_loops, dico_interactions

# fonction qui associe le code une lettre a un residu en code 3 lettres
def code_aa(dico_aa, residue):
	if dico_aa.has_key(residue):
		return dico_aa[residue]
	else:
		return "X"

	
	
def tef_version1(nom_pdb, nom_out, PDBlabel, seuil, long_min, long_max, renumber, nom_mod="", nom_chain=""):
	
	#fic_out = open(nom_out+"/"+nom_pdb.split('/')[-1].split('.')[0]+"_summary.txt", "w")
	fic_out = open(nom_out+"/"+PDBlabel+"_solutions.tef", "w")
	
	fic_pdb = open(nom_pdb, "r")
	lignes_pdb = fic_pdb.readlines()
	fic_pdb.close()
	
	# on recupere les donnees dans le fichier pdb et on calcule les distances
	# on a deja fixe le seuil min de distance et une longueur min aussi
	res_coord = coord_pdb(lignes_pdb)
	dico_coord = res_coord[0]
	dico_seqres = res_coord[1]
	dico_contigu = res_coord[2]
	dico_missing_res = res_coord[3]
	
	
	# on trie les clefs du dico pour afficher les modeles dans l'ordre
	clefs_mod = dico_coord.keys()
	clefs_mod.sort()
	
	# on effectue l'operation pour chaque modele de la structure
	for clef_mod in clefs_mod:
		# on regarde toutes les chaines aussi
		clefs_chaine = dico_coord[clef_mod].keys()
		clefs_chaine.sort()
	
		for clef_chaine in clefs_chaine:
			liste_finale = []
	
	
			if len(dico_coord) > 1:
				fic_out.write(" Modele "+clef_mod.strip())
			else:
				fic_out.write(" Modele X ")
	
			if len(dico_coord[clef_mod]) > 1:
				fic_out.write(" Chaine "+clef_chaine)
			else:
				fic_out.write(" Chaine X ")
	
			fic_out.write("\n")
	
			# on ouvre le fichier de resultat au format texte
			fic_res = open(nom_out+"/"+PDBlabel+"_model_"+clef_mod.strip()+"positions_chain_"+clef_chaine.strip()+".tef", "w")
	
			# on determine la sequence
			sequence = ""
			# si le fchier pdb est normal on a le champ SEQRES et donc tout va bien
			if dico_seqres != {}:
				seq_chaine = string.split(dico_seqres[clef_chaine])
				for res in range(len(seq_chaine)):
					sequence += code_aa(dico_aa, seq_chaine[res])
	
			# cas ou il n'y a pas de champ SEQRES dans le fichier pdb (fichier pdb genere par un logiciel de visualisation par exemple)
			else:
				for res in range(len(dico_coord[clef_mod][clef_chaine])):
					sequence += code_aa(dico_aa, dico_coord[clef_mod][clef_chaine][res][4])
	
			# on va verifier si on a le meme nombre de residus dans SEQRES et dans les coordonnees
			# si ce n'est pas le cas on force la renumerotation pour eviter certains problemes de trous
			#if dico_seqres != {} and len(dico_coord[clef_mod][clef_chaine]) >= len(sequence) and dico_contigu[clef_mod][clef_chaine] == False and not renumber:
			if len(dico_coord[clef_mod][clef_chaine]) >= len(sequence) and dico_contigu[clef_mod][clef_chaine] == False and not renumber:
				renumber = True
	
			#print str(dico_contigu[clef_mod][clef_chaine])+" %d %d<br>" %(len(dico_coord[clef_mod][clef_chaine]), len(sequence))
	
			res_calcul_ext = calcul_dist_ext(dico_coord[clef_mod][clef_chaine], seuil, long_min, long_max, renumber)
			liste_closed_loops = res_calcul_ext[0]
			dico_interactions = res_calcul_ext[1]
	
			if liste_closed_loops == []:
				fic_out.write("\n")
				fic_out.write("SEQ  "+sequence+"\n")
				fic_out.write("TEF  "+"-"*len(sequence)+"\n")
				fic_out.write("TEF  "+"-"*len(sequence)+"\n")
				fic_out.write("\n")
				continue
	
			# cas du premier fragment
			prec = [liste_closed_loops[0][0], liste_closed_loops[0][1], liste_closed_loops[0][2]]
	
			# on parcourt la liste des closed loops
			for li in range(1, len(liste_closed_loops)):
				#print "prec : "+str(prec)
				#print "liste_closed_loops[li] : "+str(liste_closed_loops[li])
	
				# on peut faire la comparaison entre les differents fragments
				# puisque prec a deja ete defini
				if liste_closed_loops[li][0]-prec[1] > -5 and liste_closed_loops[li][0] > prec[0] and liste_closed_loops[li][1] > prec[1]:
					#print prec
					#print liste_closed_loops[li][0:3]
					#print "je garde : "+str(prec)
					liste_finale.append(prec)
					prec = [liste_closed_loops[li][0], liste_closed_loops[li][1], liste_closed_loops[li][2]]
	
				else:
					if liste_closed_loops[li][2] < prec[2]:
						#print "-> meilleur"
						prec = [liste_closed_loops[li][0], liste_closed_loops[li][1], liste_closed_loops[li][2]]
	
			# on recupere le dernier tef
			liste_finale.append(prec)
	
	
			# on verifie s'il n'y a pas de trou dans la structure empechant d'assigner les TEF
			if dico_missing_res[clef_mod][clef_chaine]:
				fic_out.write("\n")
				fic_out.write("SEQ  "+sequence+"\n")
				fic_out.write("TEF  "+"-"*len(sequence)+"\n")
				fic_out.write("TEF  "+"-"*len(sequence)+"\n")
				fic_out.write("\n")
				continue
	
	
			# on rappelle que l'on a force la renumerotation a cause des differences de taille de sequence
			#if dico_seqres != {} and len(dico_coord[clef_mod][clef_chaine]) >= len(sequence) and dico_contigu[clef_mod][clef_chaine] == False and renumber == True:
	#                 if renumber == True:
	#                         fic_out.write(" Renumber T\n")
	#                 else:
	#                         fic_out.write(" Renumber F\n")
			long_moy=0.0
			for i in range(len(liste_closed_loops)):
				res_i = liste_closed_loops[i][3]+str(liste_closed_loops[i][0])
				res_j = liste_closed_loops[i][4]+str(liste_closed_loops[i][1])
				fic_res.write("%4s    %4s  %1.2lfA   %2d\n" %(res_i, res_j, liste_closed_loops[i][2], (liste_closed_loops[i][1]-liste_closed_loops[i][0]+1)))
				long_moy+=(liste_closed_loops[i][1]-liste_closed_loops[i][0]+1)
	
			long_moy/=float(len(liste_closed_loops))
	
	
			# on affiche la sequence de la structure
			fic_out.write("SEQ  "+sequence+"\n")
			fic_res.write("\n\n")
			fic_res.write(sequence+"\n")
	
			# on affiche la liste finale
			tef_haut="-"*len(sequence)
			tef_bas="-"*len(sequence)
	
			# il faut pouvoir corriger si les indices du fichier pdb ne debutent pas a 1
			if not renumber or renumber == False:
				diff_ind = dico_coord[clef_mod][clef_chaine][0][3]-1
			else:
				diff_ind = 0
			for i in range(len(liste_finale)):
				if i%2 == 0:
					tef_haut = tef_haut[0:liste_finale[i][0]-1-diff_ind] + "T"*(liste_finale[i][1]-liste_finale[i][0]+1) + tef_haut[liste_finale[i][1]-diff_ind:]
				else:
					tef_bas = tef_bas[0:liste_finale[i][0]-1-diff_ind] + "T"*(liste_finale[i][1]-liste_finale[i][0]+1) + tef_bas[liste_finale[i][1]-diff_ind:]
	
			fic_out.write("TEF  "+tef_haut+"\n")
			fic_res.write(tef_haut+"\n")
			fic_out.write("TEF  "+tef_bas+"\n")
			fic_res.write(tef_bas+"\n")
	
			# on ferme le fichier de resultat
			fic_res.close()
	
	# on ferme le fichier global de res		
	fic_out.close()




if __name__ == "__main__":
	# on verifie le bon nombre d'arguments
	if len(sys.argv) < 5:
		print "Usage: "+sys.argv[0]+" -pdb fic_pdb -o fic_out"
		print "             options : -d [float]  distance max entre les extremites"
		print "                       -min  [int]    longueur minimum d'un TEF"
		print "                       -max  [int]    longueur maximum d'un TEF"
		print "                       -renum [F ou T] si on renumerote ou pas les residus"
		print "                       -mod   [int]    numero du modele a afficher"
		print "                       -chain [char]   nom de la chaine a afficher"
		sys.exit()

	seuil = 10.0
	long_min = 10
	long_max = 100
	renumber = False
	nom_out = "./"
	
	#recuperation des arguments rentres par l'utilisateur
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-pdb":
			nom_pdb = sys.argv[i+1]
		elif sys.argv[i] == "-o":
			nom_out = sys.argv[i+1]
		elif sys.argv[i] == "-seuil":
			seuil = float(sys.argv[i+1])
		elif sys.argv[i] == "-lmin":
			long_min = int(sys.argv[i+1])
		elif sys.argv[i] == "-lmax":
			long_max = int(sys.argv[i+1])
		elif sys.argv[i] == "-renum" and string.upper(sys.argv[i+1]) == "T":
			renumber = True
		elif sys.argv[i] == "-renum" and string.upper(sys.argv[i+1]) == "F":
			renumber = False
		elif sys.argv[i] == "-mod":
			nom_mod = sys.argv[i+1]
		elif sys.argv[i] == "-chain":
			nom_chain = sys.argv[i+1]

	PDBlabel = nom_pdb.split('/')[-1].split('.')[0]

	tef_version1(nom_pdb, nom_out, PDBlabel, seuil, long_min, long_max, renumber, nom_mod, nom_chain)


