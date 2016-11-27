#!/usr/bin/python2.7
import string
from math import *

def parsePDB(infile) :
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
        
        # parcourt le PDB et cree un dictionnaire que l'on appelle dPDB
        for line in lines :
                if line[0:5] == "MODEL" :
                        conf=int(line[6:14])
                        # si la cle n'existe pas deja, on l'ajoute au dictionnaire 
                        if not conf in dPDB["conformation"] :
                                dPDB["conformation"].append(conf)
                                dPDB[conf] = {}
                                dPDB[conf]["liste_residus"] = []
                elif line[0:4] == "ATOM" :
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

def writePDB(dPDB, filout = "out_pdb.pdb") :
    """Ecrit les coordonnes de chaque atome dans un fichier .pdb"""
    # cette fonction permet d'ecrire le dictionnaire obtenu avec le parseur dans un fichier de sortie nomme out_pdb

    fout = open(filout, "w")

    for conf in dPDB["conformation"]:
        fout.write("CEMA  %4s                          %8.3f%8.3f%8.3f\n"%(conf, dPDB[conf]["XCM"], dPDB[conf]["YCM"], dPDB[conf]["ZCM"]))
        fout.write("RAGI  %4s                                                   %8.3f\n"%(conf, dPDB[conf]["rayon_giration"]))
        for res in dPDB[conf]["liste_residus"] :
            fout.write("CMRE  %4s                  %3s     %8.3f%8.3f%8.3f\n"%(conf, res, dPDB[conf][res]["XCM"], dPDB[conf][res]["YCM"], dPDB[conf][res]["ZCM"]))
            fout.write("D_CM  %4s                  %3s                              %8.3f\n"%(conf, res, dPDB[conf][res]["dist_CM"]))
            fout.write("RMSD  %4s                  %3s                              %8.3f\n"%(conf, res, dPDB[conf][res]["rmsd"]))
            for atom in dPDB[conf][res]["liste_atomes"] :
                                        fout.write("ATOM  %4s       %4s      %4s     %8.3f%8.3f%8.3f\n"%(conf, atom, res,dPDB[conf][res][atom]["x"], dPDB[conf][res][atom]["y"],dPDB[conf][res][atom]["z"]))
                    #                 			   "ATOM" CONF      atome    RES      X     Y     Z                
    fout.close()