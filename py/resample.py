#!/usr/bin/env python

"""
Tools for resampling spectra and diagonalizing the resulting covariance.

Stephen Bailey, LBL
May 2014
"""

import sys
import os
import os.path
import numpy as N
import pylab as P
from scipy import sparse
import scipy.sparse.linalg
from scipy.signal import bspline
from scipy.interpolate import LSQUnivariateSpline
from scipy.interpolate import fitpack as spl

from glob import glob
import scipy.sparse
import scipy.sparse.linalg
from scipy.sparse.construct import spdiags
from time import time, sleep

def sym_sqrt(A):
    """
    Returns Q such that A = QQ in a matrix multiplication sense,
    where A is an input Hermitian numpy ndarray.
    
    Adapted from A. Bolton, Utah
    """
    w, V = N.linalg.eigh(A)    
    wx = spdiags(N.sqrt(w), 0, len(w), len(w))
    Q = V.dot( wx.dot(V.T) )
    return Q
    
def resolution_from_icov(iCov):
    """
    Given inverse covariance matrix iCov as numpy ndarray,
    returns (R, ivar) such that
        Cov = Inverse(iCov)
        xCov = R Cov R.T = diag(1/ivar) = Diagonalized Covariance
        xiCov = Inverse(xCov) = diag(ivar)
        iCov = R.T xiCov R
        
    R is the "resolution matrix" which transforms a vector with inverse
    covariance iCov into a basis with diagonal errors given by inverse
    variances ivar.
        
    Adapted from A. Bolton, Utah
    """
    Q = sym_sqrt(iCov)
    norm = N.sum(Q, axis=1)
    R = N.outer(1.0/norm, N.ones(norm.size)) * Q
    return R, norm**2

def timeit(label, t0):
    ### print label, time() - t0
    return time()

def spline_fit(x, y, weights, xfit):
    """
    Fit a cubic b-spline to y vs. x with weights and full error propagation.
    Return fit with diagonalized errors.

    Inputs:
        x, y, weights : numpy arrays sorted by ascending x
    
    Returns 1D arrays:
        coadd_y[]
        coadd_ivar[]
        
    Where coadd_ivar is the diagonal inverse variance of coadd_y.
    """
    #- Setup cubic b-spline basis matrix: y = A coeff
    t0 = time()
    ny = len(y)
    dxfit = N.gradient(xfit)
    A = N.zeros( (ny, len(xfit)) )
    for i in range(A.shape[0]):
        A[i] = bspline((x[i] - xfit)/dxfit, 3)

    t0 = timeit("Setup A", t0)

    #- Convert to sparse matrices
    A = scipy.sparse.dia_matrix(A)
    W = spdiags(weights, 0, ny, ny)
    t0 = timeit("Sparsify A", t0)

    #- Generate inverse covariance matrix
    iCovCoeff = A.T.dot(W.dot(A))

    #- Solve "y = A coeff" with weights -> (A.T W) y = (A.T W A) coeff
    wy = A.T.dot(W.dot(y))
    t0 = timeit("Apply weights", t0)
    # coeff = scipy.sparse.linalg.lsqr(iCovCoeff, wy)[0]   
    # t0 = timeit("Solve", t0)

    #- Solve with banded matrices
    iCovDiag = iCovCoeff.todia().data[-1::-1, :]  #- have to flip for solver
    ndiag = iCovDiag.shape[0]/2
    coeff = scipy.linalg.solve_banded((ndiag, ndiag), iCovDiag, wy)
    t0 = timeit("Solve", t0)

    #- Evaluate b-spline solution at xfit
    B = N.zeros( (len(xfit), len(xfit)) )
    for i in range(B.shape[0]):
        B[i] = bspline((xfit[i] - xfit)/dxfit, 3)

    coadd = B.dot(coeff)
    t0 = timeit("Evaluate", t0)
    
    ### This is the expensive step as matrices get larger
    ### could use eig_banded and invert with that
    #- Convert iCov for coefficients into iCov for the coadd itself
    #-   coadd = B coeff
    #-   CovCoadd = B CovCoeff B.T
    #-   iCovCoadd = iB.T iCovCoeff iB,  since B is square and invertable
    iB = N.linalg.inv(B)
    iCovCoadd = iB.T.dot( iCovCoeff.dot(iB) )
    t0 = timeit("Get iCov(coadd)", t0)
    
    ### This is the next most expensive step as matrices get larger
    #- Diagonalize errors
    try:
        R, coivar = resolution_from_icov(iCovCoadd)
        Rc = R.dot(coadd)
    except:
        print "ERROR: Unable to find R for range", xfit[0], xfit[1]
        return N.zeros(len(coadd)), N.zeros(len(coadd))
    t0 = timeit("Diagonalize", t0)
        
    return Rc, coivar

