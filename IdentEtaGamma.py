#!/usr/bin/env python

import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
import math
import argparse
import cPickle
import logging

from operator import mul, div, eq, ne, add, ge, le, itemgetter
from itertools import izip
from itertools import compress
from numpy import array, log, exp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
from numpy.random import RandomState
from scipy.stats import chi2
from collections import defaultdict


def compSND(tau1,tau2):
    G1 = tau1.shape[1]
    G2 = tau2.shape[1]
        
    snd = np.zeros((G1,G2),dtype=np.int)
    N = tau1.shape[0]
    for g in range(G1):
        #snd[g,g] = 0
        for h in range(G2):
            overlap = 0.0;
            for v in range(N):
                idg = np.argmax(tau1[v,g,:])
                idh = np.argmax(tau2[v,h,:])
                if(idg == idh):
                    overlap += 1 
                
            snd[g,h] = N - overlap
                
    return snd

def main(argv):
    parser = argparse.ArgumentParser()
    
    parser.add_argument("prefix")    

    parser.add_argument("eta_file", help="gene assignments")

    parser.add_argument("var_file", help="variant frequencies")

    parser.add_argument("tau_file", help="inferred SNPs")

    parser.add_argument("gamma_file", help="relative frequency")

    args = parser.parse_args()

    eta = p.read_csv(args.eta_file, header=0, index_col=0)

    variants = p.read_csv(args.var_file, header=0, index_col=0)

    gamma = p.read_csv(args.gamma_file, header=0, index_col=0)
    gamma_matrix = gamma.as_matrix()

    tau = p.read_csv(args.tau_file, header=0, index_col=0)
    #import ipdb; ipdb.set_trace()
    tau_names = list(tau.columns.values)

    tau_names1 = tau_names[1:]

    tau_idxs = tau_names1[::4]
    tau_ints = [int(x) for x in tau_idxs]
    tau_ints_array = np.array(tau_ints)
    tau_ints_array /= 4

#,Position,4,5,6,7,12,13,14,15,20,21,22,23

    eta_matrix = eta.as_matrix()
    eta_matrix = eta_matrix[:,tau_ints_array]    
    gamma_matrix = gamma_matrix[:,tau_ints_array]
    var_matrix = variants.as_matrix()
    L = var_matrix.shape[0] 
    #import ipdb; ipdb.set_trace()
    tau_matrix = tau.as_matrix()
    tau_matrix = np.delete(tau_matrix,0,1)
    tau_matrix[tau_matrix < 0.5] = 0.0
    tau_matrix[tau_matrix >= 0.5] = 1.0

    V = tau_matrix.shape[0]
    G = tau_matrix.shape[1]/4
    
    tau_array = np.reshape(tau_matrix,(V, G,4)) 
    snd = compSND(tau_array,tau_array)
    sndP = snd/float(L)
    G = eta_matrix.shape[1]
    Z = eta_matrix.shape[0]

    for g in range(G):
    
        for h in range(g):

            etaDist = np.sum(eta_matrix[:,g] != eta_matrix[:,h])
            
            etaDist = etaDist/float(Z)

            gv = gamma_matrix[:,g]
            hv = gamma_matrix[:,h]
            gammaDist = np.dot(gamma_matrix[:,g],gamma_matrix[:,h])
    
            gammaDist /= np.sqrt(gv.dot(gv)) 
            gammaDist /= np.sqrt(hv.dot(hv))
            
            print args.prefix + "," + str(g) + "," + str(h) + "," + str(etaDist) + "," + str(sndP[g,h]) + ","+ str(gammaDist)

if __name__ == "__main__":
    main(sys.argv[1:])
