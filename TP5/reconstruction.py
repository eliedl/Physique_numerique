#!/usr/bin/env python
# -*- coding: utf-8 -*-
# TP reconstruction TDM (CT)
# Prof: Philippe Després
# programme: Dmitri Matenine (dmitri.matenine.1@ulaval.ca)


# libs
import numpy as np
import time
import matplotlib.pyplot as plt

# local files
import geometry as geo
import util as util
import CTfiltre as CTfilter
#fuckoff
def distance_point_droite(x, y, x0, y0, theta):
    dx =  (x0-x)*np.cos(theta) - (y0-y)*np.sin(theta)
    return dx

## créer l'ensemble de données d'entrée à partir des fichiers
def readInput():
    # lire les angles
    [nbprj, angles] = util.readAngles(geo.dataDir+geo.anglesFile)

    print("nbprj:",nbprj)
    print("angles min and max (rad):")
    print("["+str(np.min(angles))+", "+str(np.max(angles))+"]")

    # lire le sinogramme
    [nbprj2, nbpix2, sinogram] = util.readSinogram(geo.dataDir+geo.sinogramFile)

    if nbprj != nbprj2:
        print("angles file and sinogram file conflict, aborting!")
        exit(0)

    if geo.nbpix != nbpix2:
        print("geo description and sinogram file conflict, aborting!")
        exit(0)

    return [nbprj, angles, sinogram]


## reconstruire une image TDM en mode retroprojection
def laminogram():
    
    [nbprj, angles, sinogram] = readInput()

    # initialiser une image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox))
    # "etaler" les projections sur l'image
    # ceci sera fait de façon "voxel-driven"
    # pour chaque voxel, trouver la contribution du signal reçu
    x0 = y0 = geo.nbvox/2
    half_sino = sinogram.shape[1]/2
    for j in range(geo.nbvox): # colonnes de l'image
        #print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
                #rotation et translation de l'espace
                dist = np.array(distance_point_droite(i, j, x0, y0, angles))
                dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
                index_inf = (np.floor(half_sino + dist)).astype(int)
                index_sup = (np.ceil(half_sino + dist)).astype(int)
                residue = dist%1
                x = np.arange(0, len(sinogram))
                inf = sinogram[(x, index_inf)]
                sup = sinogram[(x, index_sup)]
                image[j,i] += np.sum(residue*(sup-inf) + inf)

    util.saveImage(image, "lam")


## reconstruire une image TDM en mode retroprojection
def backproject():
    
    [nbprj, angles, sinogram] = readInput()
    
    # initialiser une image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox))
    
    ### option filtrer ###
    #CTfilter.filterSinogram(sinogram)
    ######
    
    # "etaler" les projections sur l'image
    # ceci sera fait de façon "voxel-driven"
    # pour chaque voxel, trouver la contribution du signal reçu
    x0 = y0 = geo.nbvox / 2
    half_sino = sinogram.shape[1] / 2
    for j in range(geo.nbvox): # colonnes de l'image
        #print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
            dist = np.array(distance_point_droite(i, j, x0, y0, angles))
            dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
            index_inf = (half_sino + dist).astype(int)
            x = np.arange(0,len(sinogram))
            image[j, i] += np.sum(sinogram[(x,index_inf)])

    util.saveImage(image, "fbp")


def rotate_line(x, x0, y0, theta):
    nx = -(x - x0)*np.sin(theta) + x0
    ny = (x - x0)*np.cos(theta) + y0
    return nx, ny

## reconstruire une image TDM en mode retroprojection
def reconFourierSlice():

    [nbprj, angles, sinogram] = readInput()

    image = np.zeros((geo.nbvox, geo.nbvox), 'complex')
    fft_sino = []
    for line in sinogram:
        fft_sino.append(np.fft.fftshift(np.fft.fft(np.fft.ifftshift(line))))

    sinogram = np.array(fft_sino)

    linefunc = lambda x, a, b: a*x+b
    testat = []
    for j in range(geo.nbvox):  # colonnes de l'image
        print("working on image column: " + str(j + 1) + "/" + str(geo.nbvox))
        for i in range(geo.nbvox):

            x, y  =  -(j - geo.nbvox/2),(i - geo.nbvox/2)
            r = np.linalg.norm([x,y])#/factor
            if r == 0:
                #continue
                r = 1e-20
            theta = np.arctan2(x,y)%(2*np.pi)
            index = np.argmin(np.abs(theta- angles))
            line = sinogram[index]
            print(len(line), geo.nbvox)
            rindex = linefunc(r, -168/96, 168)
            if rindex < 0:
                continue

            rinf = (np.floor(rindex)).astype(int)
            rsup = (np.ceil(rindex)).astype(int)
            reste = rindex%1
            inf = line[rinf]
            sup = line[rsup]
            #image[j, i] += (reste * (sup - inf) + inf)
            image[j, i] += (sup+ inf)/2
            testat.append(rindex)

    image_ifft = np.fft.fftshift(np.fft.ifft2((np.fft.ifftshift(((image))))))
    util.saveImage(abs(image_ifft), "TdF slice Hi-Res")


## main ##

#laminogram()
#backproject()

reconFourierSlice()
#t = []
#for i in range(10):
#    start_time = time.time()
#    reconFourierSlice()
#    t.append(time.time() - start_time)
#t = np.array(t)
#print('Temps moyen d\'exécution', np.mean(t))



