#!/usr/bin/python2.7
import sys
from parse import *
from parse2 import *
from cenmass import *
from distances import *
from math import *
from usage import *
import matplotlib.pyplot as plt
from pylab import *
import numpy as np

# Fichier pdb a traiter : 
try:
    infile = sys.argv[sys.argv.index("-pdb")+1]                
    # infile est le premier fichier entre en argument apres l'expression -pdb
    print "pdb to treat:", infile
    
except:    
    usage() # Si le programme ne peut pas etre execute, il entre dans l'instruction except qui contient la fonction usage()
    # La fonction usage renseigne l'utilisateur du programme sur la facon dont il doit s'en servir et sur la ligne de commande a entrer pour executer le programme. 
                       
    print "ERROR: please, enter the name of the pdb input"
    sys.exit()

########################################################################################################################################################
######################################################################## MAIN ##########################################################################
########################################################################################################################################################

dico = parsePDB(infile)
centerMassOfConf(dico)
centerMassOfResidue(dico)
distance(dico)
rayon_giration(dico)
rmsd_loc(dico)
dRMSD = rmsdMoy(dico)
dMEAN = distanceMoy(dico)
writePDB(dico)
writeRayGir(dico)
writeMEAN(dMEAN)
writeRMSDMOY(dRMSD)

dico2 = parse2PDB(infile)
rmsd(dico2)
rmsd_loc(dico2)
write2PDB(dico2)
writeRMSD(dico2)

########################################################################################################################################################
###################################################################### COURBES #########################################################################
########################################################################################################################################################

## COURBE RMSD GLOBAL 

# On cree l'axe des abscisse : autant de points que de cles dans dico2 i.e. nb+1 (nb conformations a etudier + 1 conformation d'origine)
listex = []
i = 0
while i<len(dico2)-2 :
	listex.append(i)
	i += 1
x = listex[1:]

# On cree ensuite l'axe des ordonnees : les valeurs correspodant la cle rmsd pour chaque conformation
liste_RMSD = []
for conf in dico2["conformation"] :
	liste_RMSD.append(dico2[conf]["rmsd"])
rmsd = liste_RMSD[1:]
plt.figure(1)
plt.subplot(211)
plt.plot(x, rmsd)
plt.title("Courbe RMSD")
plt.xlabel("Temps")
plt.ylabel("Distance en Angstrom")

## COURBE RAYON DE GIRATION

# L'axe des abscisses reste le meme que pour la courbe du RMSD global

# On cree ensuite l'axe des ordonnees : les valeurs correspondant a la cle rayon_giration pour chaque conformation 
listeRayon = []
for conf in dico["conformation"]:
	listeRayon.append(dico[conf]["rayon_giration"])
Rayon = listeRayon[1:]
plt.subplot(212)
plt.plot(x, Rayon,"g")
ablineValues = []
for i in x :
    ablineValues.append(dico[0]["rayon_giration"])
plt.plot(x, ablineValues, "b")
plt.title("Courbe Rayon giration")
plt.xlabel("Temps")
plt.ylabel("Distance en Angstrom")

    
## COURBE DISTANCE MOYENNE AU CENTRE DE MASSE 

x_list = []
i = 0
while i < len(dMEAN) -1 : 
    x_list.append(i)
    i += 1
x_nb = x_list[1:]

mean_list = []
reslist = dMEAN["liste_residus"]
for res in reslist : 
    mean_list.append(dMEAN[res])
mean = mean_list[1:]
plt.figure(2)
plt.subplot(211)
plt.plot(x_nb, mean, "r", label = "Distance moyenne au centre de masse")
plt.title("Courbe distance moyenne au centre de masse")
plt.xlabel("Numero du residu")
plt.ylabel("Distance en Angstrom")

## COURBE RMSD MOYEN POUR CHAQUE ACIDE AMINE 

rmsd_list = []
reslist = dRMSD["liste_residus"]
for res in reslist :
    rmsd_list.append(dRMSD[res])
rmsd = rmsd_list[1:]
plt.subplot(212)
plt.plot(x_nb, rmsd, "b", label="RMSD moyen")
plt.title("Courbe RMSD moyen")
plt.xlabel("Numero du residu")
plt.ylabel("Distance en Angstrom")


## COURBES LOCALES POUR LES RESIDUS 39 ET 76

liste_RMSD_39 = []
for conf in dico["conformation"] :
	liste_RMSD_39.append(dico2[conf]["39"]["rmsd"])
rmsd_39 = liste_RMSD_39[1:]
plt.figure(3)
plt.subplot(211)
plt.plot(x, rmsd_39, "y")
plt.title("Courbe RMSD local pour l'Asparagine 39")
plt.xlabel("Temps")
plt.ylabel("Distance en Angstrom")

liste_dist_39 = []
for conf in dico["conformation"]:
	liste_dist_39.append(dico[conf]["39"]["dist_CM"])
Dist_39 = liste_dist_39[1:]
plt.subplot(212)
plt.plot(x, Dist_39,"g")
plt.title("Courbe distance au centre de masse pour l'Asparagine 39")
plt.xlabel("Temps")
plt.ylabel("Distance en Angstrom")

#plt.close(3)

liste_RMSD_76 = []
for conf in dico["conformation"] :
	liste_RMSD_76.append(dico2[conf]["33"]["rmsd"])
rmsd_76 = liste_RMSD_76[1:]
plt.figure(4)
plt.subplot(211)
plt.plot(x, rmsd_76, "y")
plt.title("Courbe RMSD local pour l'Asparagine 33")
plt.xlabel("Temps")
plt.ylabel("Distance en Angstrom")

liste_dist_76 = []
for conf in dico["conformation"]:
	liste_dist_76.append(dico[conf]["33"]["dist_CM"])
dist_76 = liste_dist_76[1:]
plt.subplot(212)
plt.plot(x, dist_76,"g")
plt.title("Courbe distance au centre de masse pour l'Asparagine 33")
plt.xlabel("Temps")
plt.ylabel("Distance en Angstrom")

#plt.close(4)

plt.show()















 