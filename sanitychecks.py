# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 12:56:19 2019

@author: prash
"""

# ==================================================================================================
# Testing basic properties of a metabolic model (aka sanity checks)
#
# Built upon COBRA toolbox sanity checks by Thiele et.al
# 
#  In this tutorial, we show how test for basic modeling properties of a metabolic 
#  model. The tutorial was developed during the construction of the generic human 
#  metabolic model, Recon 3D and can be applied to Recon 3D derived condition- 
#  and cell-type specific models, to the previous generic human reconstruction, 
#  Recon2, as well as other metabolic models.
#
# Content:
#  The tests include:
#  * leak test
#  * production of protons from nothing as well as from water, and/or oxygen alone
#  * production of matter when atp hydrolysis reaction is allowed to work but all uptakes are closed
#  * production of too much ATP from glucose under aerobic condition
#  * duplicated reactions
#  * empty colmuns in the model.rxnGeneMat
#  * the single gene deletion analysis runs smoothly
#  * ATP yield from different carbon sources
#  * metabolic objective functions
#  * flux consistency
#  * demand reactions with negative lower bound (should not occur based on definition of demand 
#    reactions)
#  * consistency of model.rev, which defines reaction reversibility, and the set values for the 
#    lower bounds on reactions. 
#  
# ==================================================================================================

#---------------------------------------------------------------------------------------------------
#load packages
#---------------------------------------------------------------------------------------------------
import cobra
#import pandas as pd

#---------------------------------------------------------------------------------------------------
#Define functions
#---------------------------------------------------------------------------------------------------
def myprogressbar(maxvalue):
    import progressbar
    bar = progressbar.ProgressBar(maxval=maxvalue, \
                widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    return bar;


#---------------------------------------------------------------------------------------------------
#Loading Models
#---------------------------------------------------------------------------------------------------
mymodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/mitocore_updated_15_04.xml")
oldmodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/mitocore_updated_16_01.xml")
biggmodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/Recon3D_bigg.xml")
reconmodel = cobra.io.load_matlab_model("D:/My Ph.D/metabolic model of mitochondria/models/Recon3DModel_Dec2017.mat")

#---------------------------------------------------------------------------------------------------
#Remove internal and unrealistic energy generating cycles
#---------------------------------------------------------------------------------------------------

#from cobra.flux_analysis.loopless import loopless_solution, add_loopless


#Check if the model leaks
mymodel_reaction_id = [rec.id for rec in mymodel.reactions]

exchanges_1 = [rec for rec in mymodel_reaction_id if rec.startswith('EX_')]
exchanges_2 = [rec for rec in mymodel_reaction_id if rec.startswith('DM_')]
exchanges_3 = [rec for rec in mymodel_reaction_id if rec.startswith('sink_')]
exchanges_4 = [rec.id for rec in mymodel.exchanges]
exchanges = exchanges_1 + exchanges_2 + exchanges_3 + exchanges_4
exchanges = list(set(exchanges))


for i in mymodel.exchanges:
    i.lower_bound = 0.0

leaky_reactions = []
for i in exchanges:
    if mymodel.reactions.get_by_id(i).flux > 1e-06:
        leaky_reactions.append(i)

if len(leaky_reactions) > 0:
    print("model leaks")
else:
    print("model leak free")
    

#Test if the model produces energy from water!
#Test for two reactions H2Ot and EX_h2o_c
#Just allow uptake of water keeping other exchange reactions closed

mymodel.reactions.get_by_id("H2Ot").lower_bound = -1000
mymodel.reactions.get_by_id("EX_h2o_c").lower_bound = -1000

atp_flux = mymodel.optimize().objective_value

if abs(atp_flux) > 1e-06:
    print('model produces energy from water!')
else:
    print("model DOES NOT produce energy from water!")


#Test if the model produces energy from water and oxygen!
#Test for H2ot, EX_h20_c, O2t and EX_o2_c

mymodel.reactions.get_by_id("H2Ot").lower_bound = -1000
mymodel.reactions.get_by_id("EX_h2o_c").lower_bound = -1000
mymodel.reactions.get_by_id("O2t").lower_bound = -1000
mymodel.reactions.get_by_id("EX_o2_c").lower_bound = -1000

atp_flux = mymodel.optimize().objective_value

if abs(atp_flux) > 1e-06:
    print('model produces energy from water and oxygen!')
else:
    print("model DOES NOT produce energy from water and oxygen!")

#Test if the model produces matter when atp demand is reversed!

mymodel.reactions.get_by_id("OF_ATP_MitoCore").lower_bound = -1000

atp_flux = mymodel.optimize().objective_value

if abs(atp_flux) > 1e-06:
    print('model produces matter when atp demand is reversed!')
else:
    print("model DOES NOT produce matter when atp demand is reversed!")
    

