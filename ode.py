#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###
# Name: Amelia & Gwyneth
# Student ID: 2289652
# Email: roseto@chapman.edu
# Course: PHYS220/MATH220/CPSC220 Fall 2018
# Assignment: CW 11
###

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
I = np.array([[0,1],[-1,0]])

def euler_1(initP, change):
    slope = I @ initP
    return initP + (change * slope)

def euler_2(N, u):
    xRange = np.arange(0, 10*np.pi, 2*np.pi/N)
    eulerApprox = np.zeros((len(xRange)+1, 2))
    change = xRange[1] - xRange[0]
    eulerApprox[0] = u
    n = 0
    for x in xRange:
        n += 1
        eulerApprox[n] = euler_1(eulerApprox[n-1], change)
    return eulerApprox

def Heuns_1(initP, change):
    Approx_2 = euler_1(initP, change)
    return initP + ((change/2)*(I @ (initP + Approx_2)))

def Heuns_2(N, u):
    xRange = np.arange(0, 10*np.pi, 2*np.pi/N)
    change = xRange[1] - xRange[0]
    heunApprox = np.zeros((len(xRange)+1, 2))
    heunApprox[0] = u
    n = 0
    for x in xRange:
        n += 1
        heunApprox[n,:] = Heuns_1(heunApprox[n-1,:], change)
    return heunApprox

def rungeKuttaSecond_1(initP, change):
    k1 = change*(I @ initP)
    k2 = change*(I @ (initP + (k1 / 2)))
    return initP + k2

def rungeKuttaSecond_2(N, u):
    xRange = np.arange(0, 10*np.pi, 2*np.pi/N)
    change = xRange[1] - xRange[0]
    RKSApprox = np.zeros((len(xRange)+1, 2))
    RKSApprox[0] = u
    n = 0
    for x in xRange:
        n+=1
        RKSApprox[n,:] = rungeKuttaSecond_1(RKSApprox[n-1,:], change)
    return RKSApprox

def rungeKuttaFourth_1(initP, change):
    k1 = change*(I @ initP)
    k2 = change*(I @ (initP + (k1 / 2)))
    k3 = change*(I @ (initP + (k2 / 2)))
    k4 = change*(I @ (initP + k3))
    return initP + (k1 + 2*k2 + 2*k3 + k4)/6

def rungeKuttaFourth_2(N, u):
    xRange = np.arange(0, 10*np.pi, 2*np.pi/N)
    change = xRange[1] - xRange[0]
    RKSApprox = np.zeros((len(xRange)+1, 2))
    RKSApprox[0] = u
    n = 0
    for x in xRange:
        n+=1
        RKSApprox[n,:] = rungeKuttaFourth_1(RKSApprox[n-1,:], change)
    return RKSApprox