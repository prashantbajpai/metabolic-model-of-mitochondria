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
#import copy

#---------------------------------------------------------------------------------------------------
#Define functions
#---------------------------------------------------------------------------------------------------
def myprogressbar(maxvalue):
    import progressbar
    bar = progressbar.ProgressBar(maxval=maxvalue, \
                widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    return bar;

def write_list_to_csv(yourlist, filename, header = ''):
    import csv
    with open(filename, 'w', newline = '') as out:
        csv_out = csv.writer(out)
        csv_out.writerow([header])
        for rows in yourlist:
            csv_out.writerow([rows])

def write_zip_to_csv(yourzip, filename, header_as_list = []):
    import csv
    with open(filename, 'w', newline = '') as out:
        csv_out = csv.writer(out)
        csv_out.writerow(header_as_list)
        for rows in yourzip:
            csv_out.writerow(rows)

def find_active_reactions(metabolite_id, model = mymodel):
    met_reactions = model.metabolites.get_by_id(metabolite_id).reactions
    met_active_reactions = [rec for rec in met_reactions if rec.flux > 1]
    return met_active_reactions;
#---------------------------------------------------------------------------------------------------
#Loading Models
#---------------------------------------------------------------------------------------------------
mymodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/mitocore_updated_15_04.xml")
oldmodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/mitocore_updated_16_01.xml")
biggmodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/Recon3D_bigg.xml")
reconmodel = cobra.io.load_matlab_model("D:/My Ph.D/metabolic model of mitochondria/models/Recon3DModel_Dec2017.mat")

mymodel.solver = 'gurobi'
mymodel.optimize()
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


for i in exchanges:
    mymodel.reactions.get_by_id(i).lower_bound = 0.0

mymodel.optimize()

leaky_reactions = []
for i in exchanges:
    if mymodel.reactions.get_by_id(i).flux > 1e-06:
        leaky_reactions.append(i)

if len(leaky_reactions) > 0:
    print("model leaks")
else:
    print("model leak free")

#Remove reactions that originally leaked in mitocore model
leaky_reactions_accepted =  ['Hct_MitoCore',
 'EX_acald_c',
 'LEUtec',
 'HDCAtr',
 'VALtec',
 'O2t',
 'HIStiDF',
 'CYStec',
 'LYStiDF',
 'ETOHt',
 'r2526',
 'r2532',
 'GLCt1r',
 'ILEtec',
 'SO3t_MitoCore',
 'r2534',
 'HCO3t_MitoCore',
 'ARGtiDF']
    
leaky_reactions = list(set(leaky_reactions) - set(leaky_reactions_accepted))

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

#===================================================================================================
#check overlapping exchange reaction and remove them

#memote report snapshot "D:/My Ph.D/metabolic model of mitochondria/mitocore_updated_15_04.xml"


#===================================================================================================
#Increase lowerbound of each exchange reaction and see which ones are affecting ATP most

imp_exchanges = []
for i in exchanges:
    mymodel.reactions.get_by_id(i).lower_bound = -1000.0
    atp_flux = mymodel.optimize().objective_value
    if atp_flux > 1e-06:
        imp_exchanges.append([i, atp_flux])
    mymodel.reactions.get_by_id(i).lower_bound = 0.0
    
#analyse bounds of reaction OF_ATP_MitoCore atp_c + h2o_c --> adp_c + biomass_c + h_c + pi_c 
    
mymodel.reactions.Biomass_MitoCore.bounds = -1000.0, 0.0
mymodel.optimize()   #ATP flux is thousand even though all other exchange reactions are closed

#pi_c flux
pi_active_reactions = find_active_reactions('pi_c')
#PYNP2r pi_c + uri_c <=> r1p_c + ura_c (-1000.0, 1000.0) 1000.0

#h_c
hc_active_reactions = find_active_reactions('h_c')
#EX_h_c h_c -->  (0.0, 1000.0) 1000.0
#OF_ATP_MitoCore atp_c + h2o_c <=> adp_c + biomass_c + h_c + pi_c (-1000, 1000.0) 1000.0
#r1088  <=> cit_c + h_c (-1000.0, 1000.0) 1000.0
#PRPPS atp_c + r5p_c <=> amp_c + h_c + prpp_c (-1000.0, 1000.0) 1000.0

#adp_c
adpc_active_reactions = find_active_reactions('adp_c')
#OF_ATP_MitoCore atp_c + h2o_c <=> adp_c + biomass_c + h_c + pi_c (-1000, 1000.0) 1000.0

#h2o_c
h2oc_active_reactions = find_active_reactions('h2o_c')
#EX_h2o_c h2o_c -->  (0.0, 1000.0) 1000.0
#OF_ATP_MitoCore atp_c + h2o_c <=> adp_c + biomass_c + h_c + pi_c (-1000, 1000.0) 1000.0
#H2Ot  --> h2o_c (0.0, 1000.0) 999.981

#atp_c
atpc_active_reactions = find_active_reactions('atp_c')
#OF_ATP_MitoCore atp_c + h2o_c <=> adp_c + biomass_c + h_c + pi_c (-1000, 1000.0) 1000.0
#PRPPS atp_c + r5p_c <=> amp_c + h_c + prpp_c (-1000.0, 1000.0) 1000.0

#biomass_c
biomassc_active_reactions = find_active_reactions('biomass_c')
#OF_ATP_MitoCore atp_c + h2o_c <=> adp_c + biomass_c + h_c + pi_c (-1000, 1000.0) 1000.0

testmodel = mymodel.copy()

