# -*- coding: utf-8 -*-
"""
Created on Mon Apr  5 09:15:56 2021
@author: Brisid
"""

import matplotlib.pyplot as plt

class Checks:
    from MKurbature import momentCurvature, epsilon, concreteStress, steelStress
    
    def kontroll(self):
        N = 1200 #kN, compression positive.
        b=300
        h=500
        c=35 #cover
        tol = 0.05 #measured in millimeters!
        maxIterations = 1000
        sec = [[0, 0], [b, 0], [b, h], [0,h]]
        secC = [[c, c], [b-c, 0+c], [b-c, h-c], [0+c,h-c]]
        steelBars = [[30,50,3.14*50*50*0.25],[20,450,3.14*50*50*0.25]]
        concreteUC = [40, 0.0022, 20, 0.0035]
        concreteC = [50, 0.0032, 35, 0.01]
        steel = [500, 0.01, 0.05]
        steps = 30
        fibersCore=30
        fibersCover=4
        
        k, M, outFS, outCORE, outCOVER, crushing, crushingCore, yielding, rupture, dropped, converged, convergedSteps = self.momentCurvature(sec, secC, steelBars, concreteUC, concreteC, steel, N, steps, tol, maxIterations, fibersCore, fibersCover)
 
        plt.plot(k, M)
        plt.show()
                
        print(f"Out of {steps} steps, {convergedSteps} converged.")
        

check = Checks()
check.kontroll()

