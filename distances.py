#!/usr/bin/python2.7
import string
from math import *

########################################################################################################################################################
################################################################# DISTANCES GLOBALES ###################################################################
########################################################################################################################################################

def rayon_giration(dPDB) :
    # La fonction rayon_giration permet de calculer le rayon de giration de chaque conformation du dictionnaire donne en entree
    # Le rayon de giration est la distance pour chaque conformation entre son centre de masse et l'atome le plus eloigne de ce centre de masse.
    
    conflist = dPDB["conformation"]
    for conf in conflist :
        raymax = 0
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist :
            atomlist = dPDB[conf][res]["liste_atomes"]
            for atom in atomlist :
                dist = sqrt((dPDB[conf][res][atom]["x"]-dPDB[conf]["XCM"])**2 +(dPDB[conf][res][atom]["y"]-dPDB[conf]["YCM"])**2 +(dPDB[conf][res][atom]["z"]-dPDB[conf]["ZCM"])**2 )
                if (dist >= raymax) :
                    raymax = dist
        dPDB[conf]["rayon_giration"] = raymax


def rmsd(dPDB) :
    # La fonction rmsd permet de calculer le RMSD global de chaque conformation par rapport a la conformation d'origine 
    # Le rmsd global est calcule en tenant compte des atomes CA, C, O et N uniquement.  
    
    conflist = dPDB["conformation"]
    for conf in conflist :
        n = 0
        somme = 0
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist:
            atomlist = dPDB[conf][res]["liste_atomes"]
            for atom in atomlist :
                dist = (dPDB[conf][res][atom]["x"]-dPDB[0][res][atom]["x"])**2 + (dPDB[conf][res][atom]["y"]-dPDB[0][res][atom]["y"])**2 + (dPDB[conf][res][atom]["z"] - dPDB[0][res][atom]["z"])**2
                somme += dist
                n += 1
        rmsd = sqrt(somme/n)
        dPDB[conf]["rmsd"] = rmsd

########################################################################################################################################################
################################################################# DISTANCES LOCALES ####################################################################
########################################################################################################################################################
        
def rmsd_loc(dPDB) :
    # La fonction rmsd_loc permet de calculer le RMSD local de chaque residu entre une conformation et la confomation d'origine 
    # Le rmsd local est calcule en tenant compte de tous les atomes pour chaque residu. 
    
    conflist = dPDB["conformation"]
    for conf in conflist :
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist:
            n = 0
            somme = 0
            atomlist = dPDB[conf][res]["liste_atomes"]
            for atom in atomlist :
                dist = (dPDB[conf][res][atom]["x"]-dPDB[0][res][atom]["x"])**2 + (dPDB[conf][res][atom]["y"]-dPDB[0][res][atom]["y"])**2 + (dPDB[conf][res][atom]["z"] - dPDB[0][res][atom]["z"])**2
                somme += dist
                n += 1
            dPDB[conf][res]["rmsd"] = sqrt(somme/n)
            
def distance(dPDB) :
    # La fonction distance permet de calculer la distance entre le centre de masse d'un residu et le centre de masse de la conformation en question 
    
    conflist = dPDB["conformation"]
    for conf in conflist :
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist :
            dist = sqrt((dPDB[conf][res]["XCM"]-dPDB[conf]["XCM"])**2 +(dPDB[conf][res]["YCM"]-dPDB[conf]["YCM"])**2 +(dPDB[conf][res]["ZCM"]-dPDB[conf]["ZCM"])**2 )
            dPDB[conf][res]["dist_CM"] = dist

########################################################################################################################################################
################################################################# DISTANCES MOYENNES ###################################################################
########################################################################################################################################################

def distanceMoy(dPDB) :
    # La fonction distanceMoy permet de calculer la distance moyenne au cours du temps entre le centre de masse d'un residu et le centre de masse de la proteine 
    
    dMEAN = {}
    dMEAN["liste_residus"] = []
    n = 0
    conflist = dPDB["conformation"]
    for conf in conflist :
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist :
            if not res in dMEAN["liste_residus"] : 
                dMEAN["liste_residus"].append(res)
                dMEAN[res] = dPDB[conf][res]["dist_CM"]
            else :
                dMEAN[res] += dPDB[conf][res]["dist_CM"]
        n += 1
    reslist = dMEAN["liste_residus"]
    for res in reslist :
        dMEAN[res] = dMEAN[res]/n
        
    return dMEAN
    
def rmsdMoy(dPDB) :
    # La fonction rmsdMoy permet de calculer le RMSD local moyen au cours du temps pour chaque acide amine 
    
    dRMSD = {}
    dRMSD["liste_residus"] = []
    n = 0
    conflist = dPDB["conformation"]
    for conf in conflist : 
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist :
            if not res in dRMSD["liste_residus"] :
                dRMSD["liste_residus"].append(res)
                dRMSD[res] = dPDB[conf][res]["rmsd"]
            else : 
                dRMSD[res] += dPDB[conf][res]["rmsd"]
        n += 1
    reslist = dRMSD["liste_residus"]
    for res in reslist :
        dRMSD[res] = dRMSD[res]/n
        
    return dRMSD
        
            
########################################################################################################################################################
############################################################# AFFICHAGE DANS FICHIERS DE SORTIE ########################################################
########################################################################################################################################################

def writeRayGir(dPDB, fileout = "giration.pdb") :
    
    fout = open(fileout, "w")
    
    conflist = dPDB["conformation"]
    for conf in conflist :
        fout.write("CONFO %4s  RAYON_GIRATION %8.3f\n"%(conf, dPDB[conf]["rayon_giration"]))
        #          "CONFO" conf  "RAYON_GIRATION" rayon
    
    fout.close()
    
def writeRMSD(dPDB, fileout = "rmsd.pdb") :
    
    fout = open(fileout, "w")
    
    conflist = dPDB["conformation"]
    for conf in conflist :
        fout.write("CONFO %4s    RMSD    %8.3f\n"%(conf, dPDB[conf]["rmsd"]))
        #          "CONFO" conf  "RMSD"  rmsd
    
    fout.close()
    
def writeMEAN(dMEAN, fileout = "dist_mean.pdb") :
    
    fout = open(fileout, "w")
    
    reslist = dMEAN["liste_residus"]
    for res in reslist :
        fout.write("RESID       %4s     Dist_moy        %8.3f\n"%(res, dMEAN[res]))
        #          "RESID"      res    "Distance_moy"  distance
        
    fout.close()
    
def writeRMSDMOY(dRMSD, fileout = "rmsd_mean.pdb") :
    
    fout = open(fileout, "w")
    
    reslist = dRMSD["liste_residus"]
    for res in reslist :
        fout.write("RESID       %4s     RMSD_Moy        %8.3f\n"%(res, dRMSD[res]))
        #          "RESID"     res     "RMSD_Moy"       rmsd
        
    fout.close()