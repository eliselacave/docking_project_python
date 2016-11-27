#!/usr/bin/python2.7
import string
from math import *

def parse2PDB(infile) :
    # Ce parseur (qui fonctionne de la meme maniere que le fichier parse.pdb) permet de parser le fichier pdb en ne retenant que les atomes du 
    # squelette peptidique de chaque residu (CA, C, O, N)
    
        # lecture du fichier PDB 
        f = open(infile, "r")
        lines = f.readlines()
        f.close()
        
        # var init
        conf= True
        firstline = True
        prevres = None
        dPDB = {}
        dPDB["conformation"] = []
        dPDB["liste_residus"] = []
        
        # parcoure le PDB
        for line in lines :
                if line[0:5] == "MODEL" :
                        conf=int(line[6:14])
                        if not conf in dPDB["conformation"] :
                                dPDB["conformation"].append(conf)
                                dPDB[conf] = {}
                                dPDB[conf]["liste_residus"] = []
                elif line[0:4] == "ATOM" :
                    if line[13:15] == "CA" or line[13:15] == "N " or line[13:15] == "C " or line[13:15] == "O " :
                        # On ne retient que les atomes du squelette peptidique
                        curres = "%s"%(line[22:26]).strip()
                        if not curres in dPDB[conf]["liste_residus"] :
                                dPDB[conf]["liste_residus"].append(curres)
                                dPDB[conf][curres] = {}
                                dPDB[conf][curres]["nom_resid"] = string.strip(line[17:20])
                                dPDB[conf][curres]["liste_atomes"] = []
                                
                        atomtype = string.strip(line[12:16])
                        dPDB[conf][curres]["liste_atomes"].append(atomtype)
                        dPDB[conf][curres][atomtype] = {}
                        dPDB[conf][curres][atomtype]["x"] = float(line[30:38])
                        dPDB[conf][curres][atomtype]["y"] = float(line[38:46])
                        dPDB[conf][curres][atomtype]["z"] = float(line[46:54])
                       
                       
                        
        return dPDB
    
def write2PDB(dPDB, filout = "out_pdb_CACON.pdb") :
    """Ecrit les coordonnees de chaque atome dans un fichier .pdb"""

    fout = open(filout, "w")

    for conf in dPDB["conformation"]:
        for res in dPDB[conf]["liste_residus"] :
            fout.write("RMSD        %8.3f\n"%(dPDB[conf][res]["rmsd"]))
            for atom in dPDB[conf][res]["liste_atomes"] :
                fout.write("ATOM  %4s      %4s     %8.3f%8.3f%8.3f\n"%(conf, res,dPDB[conf][res][atom]["x"], dPDB[conf][res][atom]["y"],dPDB[conf][res][atom]["z"]))
                    #      "ATOM" CONF     RES       X    Y   Z                
    fout.close()