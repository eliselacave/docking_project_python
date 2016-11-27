#!/usr/bin/python2.7
from math import *

########################################################################################################################################################
############################################################### POUR UNE CONFORMATION ##################################################################
########################################################################################################################################################

def centerMassOfConf(dPDB):
    """Calcule le centre de masse d'une conformation en tenant compte de tous ses residus et de tous les atomes de chaque residu"""
    # Le centre de masse de la conformation est calcule en tenant compte de tous les atomes de chaque residu.
        
    conflist = dPDB["conformation"]
        
    for conf in conflist : 
        x = y = z = length = 0.0
        
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist : 
            atlist = dPDB[conf][res]["liste_atomes"]
            for atom in atlist :
                    #at = atom[0]
                    #masseof = masse[atom]
                x += dPDB[conf][res][atom]["x"] #* masseof
                y += dPDB[conf][res][atom]["y"] #* masseof
                z += dPDB[conf][res][atom]["z"] #* masseof
                length += 1

        Xcm = float(x)/length
        Ycm = float(y)/length
        Zcm = float(z)/length
        dPDB[conf]["XCM"] = Xcm
        dPDB[conf]["YCM"] = Ycm
        dPDB[conf]["ZCM"] = Zcm
        
########################################################################################################################################################
################################################################### POUR UN RESIDU #####################################################################
########################################################################################################################################################
            
def centerMassOfResidue(dPDB):
    """Calcule le centre de masse d'un residu en tenant compte de tous ses atomes"""
        
    conflist = dPDB["conformation"]
        
    for conf in conflist : 
        
        reslist = dPDB[conf]["liste_residus"]
        for res in reslist : 
            x = y = z = length = 0.0
            atlist = dPDB[conf][res]["liste_atomes"]
            for atom in atlist :
                x += dPDB[conf][res][atom]["x"] 
                y += dPDB[conf][res][atom]["y"]
                z += dPDB[conf][res][atom]["z"] 
                length += 1
            Xcm = float(x)/length
            Ycm = float(y)/length
            Zcm = float(z)/length
            dPDB[conf][res]["XCM"] = Xcm
            dPDB[conf][res]["YCM"] = Ycm
            dPDB[conf][res]["ZCM"] = Zcm
                

