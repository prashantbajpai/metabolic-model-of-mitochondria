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

def find_active_reactions(metabolite_id, model, reaction_flux = 1):
    met_reactions = model.metabolites.get_by_id(metabolite_id).reactions
    met_active_reactions = [rec for rec in met_reactions if rec.flux > reaction_flux]
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
testmodel.optimize()

def met_critcal_reactions(met, model):
    all_reactions = find_active_reactions(met, model=model, reaction_flux=200)
    for i in all_reactions:
        import re
        x = re.findall('EX|Sink|DM', i.id)
        if len(x)==1:
            print('bounds of', i.id, 'changed from', i.bounds, 'to 0.0, 0.0')
            i.bounds = 0.0, 0.0
            model.optimize()
            atp_flux = model.optimize().objective_value
            if atp_flux < 1000.0:
                print(i.id, 'is a critical reaction');break            
        elif 'OF_ATP_MitoCore' not in i.id:
            print('bounds of', i.id, 'changed from', i.bounds, 'to 0.0, 0.0')
            i.bounds = 0.0, 0.0
            model.optimize()
            atp_flux = model.optimize().objective_value
            if atp_flux < 1000.0:
                print(i.id, 'is a critical reaction');break   
            
            
from cobra.flux_analysis import flux_variability_analysis
  
testmodel_fva = flux_variability_analysis(testmodel, pfba_factor=1.1)  
critical_reaction_id = testmodel_fva.sort_values(by = 'maximum', ascending=False).index[1:30]

critical_reactions = []
for i in critical_reaction_id:
    flux = testmodel.reactions.get_by_id(i).flux
    if flux > 100:
        critical_reactions.append(i)

for i in testmodel.reactions:
    if i.id != 'OF_ATP_MitoCore':
        if i.id != 'Biomass_MitoCore':
            i.lower_bound = 0.0
            testmodel.optimize()
            atp_flux = testmodel.optimize().objective_value
            if atp_flux < 950.0:
                print(i.id, 'is a critical reaction');break   

#FACOAL240 is a critical reaction (830) when reducing bounds to 0.0
#RE2078M is a critical reaction (1530) when reducing lower bound to 0.0
testmodel2 = testmodel.copy()

for i in range(len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 50.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break

#31 , ICDHyrm is a critical reaction

#IN FLUXES         OUT FLUXES          OBJECTIVES
#----------------  ------------------  --------------------
#fol_c      1e+03  co2_c        1e+03  OF_ATP_MitoCore  141
#h2o_c    500      dhf_c        1e+03
#o2_c      19.8    h_c        986
#hdca_c     1      succ_c     467
#glc_D_c    0.9    biomass_c  141
#ser_L_c    0.017  pyr_c       23
#phe_L_c    0.012  ac_c         6.66
#tyr_L_c    0.011  nh4_c        0.02
#asn_L_c    0.01   ptrc_c       0.007
#arg_L_c    0.007  urea_c       0.007
#                  but_c        0.006
        
for i in range(len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 150.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 34 , SUCOASm is a critical reaction
# 
# testmodel2.summary()
# IN FLUXES          OUT FLUXES          OBJECTIVES
# -----------------  ------------------  --------------------
# h2o_c      1e+03   co2_c        1e+03  OF_ATP_MitoCore  192
# fol_c    715       h_c          1e+03
# fe2_c     72.5     succ_c     675
# o2_c      19.8     thf_c      646
# gln_L_c    2.31    hco3_c     348
# hdca_c     1       biomass_c  192
# glc_D_c    0.9     dhf_c       69
# lys_L_c    0.03    ac_c         8.25
# ser_L_c    0.017   nh4_c        2.4
# phe_L_c    0.012   2oxoadp_c    0.03
# thr_L_c    0.012   acald_c      0.012
# tyr_L_c    0.011   urcan_c      0.01
# val_L_c    0.011
# chol_c     0.0105
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

for i in range(len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 200.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# #38 , CII_MitoCore is a critical reaction
# #
# #testmodel2.summary()
# #IN FLUXES          OUT FLUXES          OBJECTIVES
# #-----------------  ------------------  --------------------
# #h2o_c      1e+03   co2_c        1e+03  OF_ATP_MitoCore  201
# #fol_c    579       h_c          1e+03
# #mal_L_c  229       succ_c     936
# #o2_c      19.8     thf_c      579
# #gln_L_c    2.31    hco3_c     411
# #hdca_c     1       biomass_c  201
# #glc_D_c    0.9     ac_c         8
# #lys_L_c    0.03    nh4_c        2.4
# #ser_L_c    0.017   2oxoadp_c    0.03
# #phe_L_c    0.012   acac_c       0.023
# #thr_L_c    0.012   acald_c      0.012
# #tyr_L_c    0.011   urcan_c      0.01
# #val_L_c    0.011
# #chol_c     0.0105
# #asn_L_c    0.01
# #his_L_c    0.01
# #arg_L_c    0.007
# #met_L_c    0.005
# #cys_L_c    0.001
# =============================================================================

for i in range(len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 210.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break

# =============================================================================
# 54 , FUM is a critical reaction
# 
# testmodel2.summary()
# IN FLUXES        OUT FLUXES          OBJECTIVES
# ---------------  ------------------  --------------------
# h2o_c     1e+03  akg_c        1e+03  OF_ATP_MitoCore  977
# mal_L_c   1e+03  hco3_c     989
# gln_L_c  28      h_c        987
# o2_c     19.8    biomass_c  977
# fol_c    14      succ_c     975
# hdca_c    1      co2_c       37
# glc_D_c   0.9    nh4_c       28
# lys_L_c   0.03   thf_c       14
# ser_L_c   0.017  ac_c         8
# phe_L_c   0.012  cit_c        1.06
# thr_L_c   0.012  gly_c        0.021
# tyr_L_c   0.011  acald_c      0.012
# val_L_c   0.011  HC00250_c    0.001
# asn_L_c   0.01
# his_L_c   0.01
# arg_L_c   0.007
# met_L_c   0.005
# cys_L_c   0.001
# =============================================================================

testmodel2.reactions.FUM.bounds = -100.0, 1000.0
testmodel2.optimize()

for i in range(55,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 300.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break

# =============================================================================
# 59 , ASPTA is a critical reaction
# testmodel2.summary()
# IN FLUXES         OUT FLUXES           OBJECTIVES
# ----------------  -------------------  --------------------
# h2o_c      1e+03  akg_c         1e+03  OF_ATP_MitoCore  723
# mal_L_c  621      co2_c         1e+03
# gln_L_c  520      nh4_c         1e+03
# fol_c    430      succ_c      760
# chol_c   140      h_c         752
# ps_hs_c  140      biomass_c   723
# o2_c      19.8    pyr_c       141
# hco3_c     1.71   cit_c       141
# glc_D_c    0.9    pchol_hs_c  140
# lys_L_c    0.03   2oxoadp_c     0.03
# ser_L_c    0.017  34hpp_c       0.011
# phe_L_c    0.012  urcan_c       0.01
# thr_L_c    0.012  tsul_c        0.001
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# met_L_c    0.005
# cys_L_c    0.001
# so3_c      0.001
# =============================================================================

testmodel2.reactions.get_by_id('ASPTA').bounds = -100.0, 1000.0
testmodel2.optimize()

for i in range(60,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 380.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break

# =============================================================================
# 95 , r0791 is a critical reaction
# testmodel2.summary()
# IN FLUXES         OUT FLUXES          OBJECTIVES
# ----------------  ------------------  --------------------
# h2o_c      1e+03  co2_c        1e+03  OF_ATP_MitoCore  380
# fol_c    615      h_c          1e+03
# gln_L_c  102      succ_c     630
# o2_c      19.8    biomass_c  380
# hdca_c     1      cit_c      219
# glc_D_c    0.9    nh4_c      202
# lys_L_c    0.03   hco3_c     156
# ser_L_c    0.017  ac_c         8.99
# phe_L_c    0.012  34hpp_c      0.011
# tyr_L_c    0.011  gly_c        0.009
# val_L_c    0.011  but_c        0.006
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

for i in range(60,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 390.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break

# =============================================================================
# 214 , r0074 is a critical reaction
# testmodel2.summary()
# IN FLUXES         OUT FLUXES          OBJECTIVES
# ----------------  ------------------  --------------------
# h2o_c      1e+03  co2_c        1e+03  OF_ATP_MitoCore  393
# fol_c    594      succ_c       1e+03
# ppa_c    184      h_c        989
# gln_L_c  182      biomass_c  393
# o2_c      19.8    nh4_c      282
# hdca_c     1      hco3_c     217
# glc_D_c    0.9    ac_c        73.1
# lys_L_c    0.03   34hpp_c      0.011
# ser_L_c    0.017  gly_c        0.009
# phe_L_c    0.012  but_c        0.006
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

for i in range(60,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 400.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 219 , ORNTArm is a critical reaction
# testmodel2.summary()
# IN FLUXES         OUT FLUXES          OBJECTIVES
# ----------------  ------------------  --------------------
# h2o_c      1e+03  co2_c        1e+03  OF_ATP_MitoCore  773
# so3_c      1e+03  so4_c        1e+03
# gln_L_c  753      succ_c       1e+03
# fol_c    140      biomass_c  773
# o2_c      19.8    orn_c      565
# h_c       18.2    ac_c       392
# ppa_c     15.8    nh4_c      375
# hco3_c     1.71   pyr_c        1.83
# hdca_c     1      2oxoadp_c    0.03
# glc_D_c    0.9    34hpp_c      0.011
# lys_L_c    0.03   gly_c        0.009
# ser_L_c    0.017
# phe_L_c    0.012
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

testmodel2.reactions.get_by_id('ORNTArm').bounds = -100.0, 1000.0
testmodel2.optimize()

for i in range(220,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 475.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 233 , r0178 is a critical reaction
# testmodel2.summary()
# IN FLUXES         OUT FLUXES          OBJECTIVES
# ----------------  ------------------  --------------------
# fol_c      1e+03  co2_c        1e+03  OF_ATP_MitoCore  488
# h2o_c      1e+03  succ_c       1e+03
# so3_c      1e+03  so4_c        1e+03
# ppa_c    467      h_c        830
# gln_L_c  383      dhf_c      780
# mal_L_c  218      4abut_c    551
# o2_c      19.8    biomass_c  488
# hco3_c     1.71   thf_c      220
# hdca_c     1      orn_c      100
# glc_D_c    0.9    pyr_c        2.4
# lys_L_c    0.03   2oxoadp_c    0.03
# ser_L_c    0.017  2hb_c        0.012
# phe_L_c    0.012  34hpp_c      0.011
# thr_L_c    0.012  urcan_c      0.01
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

for i in range(220,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 500.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 296 , CKc is a critical reaction
# testmodel2.summary()
# IN FLUXES          OUT FLUXES          OBJECTIVES
# -----------------  ------------------  --------------------
# fol_c      1e+03   co2_c        1e+03  OF_ATP_MitoCore  508
# h2o_c      1e+03   succ_c       1e+03
# so3_c      1e+03   so4_c        1e+03
# ppa_c    387       h_c        857
# gln_L_c  374       dhf_c      806
# mal_L_c  281       biomass_c  508
# o2_c      19.8     4abut_c    489
# hco3_c     1.71    thf_c      194
# hdca_c     1       orn_c      100
# glc_D_c    0.9     pyr_c        2.42
# lys_L_c    0.03    2oxoadp_c    0.03
# ser_L_c    0.017   acald_c      0.012
# phe_L_c    0.012   34hpp_c      0.011
# thr_L_c    0.012
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# chol_c     0.0055
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

for i in range(296,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 520.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
 
# =============================================================================
# 344 , r0913 is a critical reaction
# testmodel2.summary()
# IN FLUXES         OUT FLUXES            OBJECTIVES
# ----------------  --------------------  ----------------------
# h2o_c      1e+03  biomass_c    1e+03    OF_ATP_MitoCore  1e+03
# so3_c      1e+03  co2_c        1e+03
# fol_c    950      succ_c       1e+03
# ppa_c    366      so4_c        1e+03
# gln_L_c  231      h_c        962
# chol_c   131      thf_c      562
# ps_hs_c  123      dhf_c      388
# mal_L_c   26.6    nh4_c      131
# o2_c      19.8    pyr_c      125
# hco3_c     1.71   orn_c      100
# hdca_c     1      4abut_c     26.7
# glc_D_c    0.9    2oxoadp_c    0.03
# lys_L_c    0.03   acald_c      0.012
# ser_L_c    0.017  34hpp_c      0.011
# phe_L_c    0.012  ppbng_c      0.00183
# thr_L_c    0.012
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

testmodel2.reactions.get_by_id('r0913').bounds = -50.0, 1000.0
testmodel2.optimize()

for i in range(345,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 560.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 346 , r0917 is a critical reaction
# testmodel2.summary()
# IN FLUXES          OUT FLUXES          OBJECTIVES
# -----------------  ------------------  ----------------------
# fol_c      1e+03   biomass_c    1e+03  OF_ATP_MitoCore  1e+03
# h2o_c      1e+03   co2_c        1e+03
# so3_c      1e+03   so4_c        1e+03
# ppa_c    409       succ_c       1e+03
# mal_L_c  315       h_c        869
# gln_L_c  216       thf_c      612
# o2_c      19.8     dhf_c      388
# hco3_c     1.71    pyr_c      210
# hdca_c     1       orn_c      100
# glc_D_c    0.9     4abut_c     99.4
# lys_L_c    0.03    2oxoadp_c    0.03
# ser_L_c    0.017   acald_c      0.012
# phe_L_c    0.012   34hpp_c      0.011
# thr_L_c    0.012
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# arg_L_c    0.007
# chol_c     0.0055
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

testmodel2.reactions.get_by_id('r0917').bounds = -50.0, 1000.0
testmodel2.optimize()

for i in range(347,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 600.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 364 , r2398B_MitoCore is a critical reaction
# testmodel2.summary()
# IN FLUXES          OUT FLUXES          OBJECTIVES
# -----------------  ------------------  ----------------------
# fol_c      1e+03   biomass_c    1e+03  OF_ATP_MitoCore  1e+03
# h2o_c      1e+03   co2_c        1e+03
# so3_c      1e+03   so4_c        1e+03
# ppa_c    384       succ_c       1e+03
# gln_L_c  232       h_c        883
# mal_L_c  108       dhf_c      554
# o2_c      19.8     thf_c      446
# hco3_c     1.71    4abut_c    120
# hdca_c     1       ptrc_c     100
# glc_D_c    0.9     pyr_c        2.43
# lys_L_c    0.03    2oxoadp_c    0.03
# ser_L_c    0.017   acac_c       0.023
# phe_L_c    0.012   acald_c      0.012
# thr_L_c    0.012
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# trp_L_c    0.008
# arg_L_c    0.007
# chol_c     0.0055
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

testmodel2.reactions.get_by_id('r2398B_MitoCore').bounds = -50.0, 1000.0
testmodel2.optimize()

for i in range(365,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 650.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 365 , r2402B_MitoCore is a critical reaction
# testmodel2.summary()
# IN FLUXES          OUT FLUXES          OBJECTIVES
# -----------------  ------------------  ----------------------
# fol_c      1e+03   biomass_c    1e+03  OF_ATP_MitoCore  1e+03
# h2o_c      1e+03   co2_c        1e+03
# so3_c    846       succ_c       1e+03
# ppa_c    376       h_c        956
# gln_L_c  222       so4_c      846
# mal_L_c  150       thf_c      504
# fe2_c     74.9     dhf_c      496
# o2_c      19.8     4abut_c    216
# hco3_c     1.71    orn_c      100
# hdca_c     1       akg_c       72.9
# glc_D_c    0.9     pyr_c        2.43
# lys_L_c    0.03    2oxoadp_c    0.03
# ser_L_c    0.017   acald_c      0.012
# phe_L_c    0.012
# thr_L_c    0.012
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# trp_L_c    0.008
# chol_c     0.0055
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

testmodel2.reactions.get_by_id('r2402B_MitoCore').bounds = -50.0, 1000.0
testmodel2.optimize()

for i in range(366,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 680.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 366 , LYStmB_MitoCore is a critical reaction
# testmodel2.summary()
# IN FLUXES          OUT FLUXES           OBJECTIVES
# -----------------  -------------------  ----------------------
# fol_c      1e+03   biomass_c     1e+03  OF_ATP_MitoCore  1e+03
# h2o_c      1e+03   co2_c         1e+03
# so3_c      1e+03   h_c           1e+03
# ppa_c    350       succ_c        1e+03
# gln_L_c  196       so4_c         1e+03
# mal_L_c  150       thf_c       598
# o2_c      19.8     dhf_c       402
# hco3_c     1.71    orn_c       100
# hdca_c     1       akg_c        62.4
# glc_D_c    0.9     pyr_c         2.42
# lys_L_c    0.03    Lpipecol_c    0.03
# ser_L_c    0.017   acald_c       0.012
# phe_L_c    0.012   urcan_c       0.01
# thr_L_c    0.012
# tyr_L_c    0.011
# val_L_c    0.011
# chol_c     0.0105
# asn_L_c    0.01
# his_L_c    0.01
# trp_L_c    0.008
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================
        
testmodel2.reactions.get_by_id('LYStmB_MitoCore').bounds = -50.0, 1000.0
testmodel2.optimize()

for i in range(367,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 720.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 368 , ARGtmB_MitoCore is a critical reaction
# testmodel2.summary()
# IN FLUXES          OUT FLUXES          OBJECTIVES
# -----------------  ------------------  ----------------------
# fol_c      1e+03   biomass_c    1e+03  OF_ATP_MitoCore  1e+03
# h2o_c      1e+03   co2_c        1e+03
# so3_c      1e+03   h_c          1e+03
# gln_L_c  248       so4_c        1e+03
# ppa_c    240       succ_c     878
# mal_L_c  234       dhf_c      676
# o2_c      19.8     thf_c      324
# hco3_c     1.71    4abut_c    266
# hdca_c     1       cit_c      156
# glc_D_c    0.9     orn_c      100
# lys_L_c    0.03    pyr_c        1.85
# ser_L_c    0.017   2oxoadp_c    0.03
# phe_L_c    0.012   acac_c       0.023
# thr_L_c    0.012   acald_c      0.012
# tyr_L_c    0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# trp_L_c    0.008
# chol_c     0.0055
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

testmodel2.reactions.get_by_id('ARGtmB_MitoCore').bounds = -50.0, 1000.0
testmodel2.optimize()

for i in range(369,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 770.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 370 , PYRt2m is a critical reaction
# testmodel2.summary()
# IN FLUXES         OUT FLUXES           OBJECTIVES
# ----------------  -------------------  --------------------
# fol_c      1e+03  co2_c         1e+03  OF_ATP_MitoCore  775
# h2o_c      1e+03  h_c           1e+03
# so3_c      1e+03  succ_c        1e+03
# ppa_c    290      so4_c         1e+03
# mal_L_c  272      biomass_c   775
# gln_L_c  241      dhf_c       580
# etoh_c    51.8    thf_c       410
# o2_c      19.8    4abut_c     274
# chol_c     9.9    orn_c       100
# pe_hs_c    3.3    cit_c        82.9
# hco3_c     1.71   acald_c      51.8
# hdca_c     1      pyr_c        21.6
# glc_D_c    0.9    mlthf_c       9.88
# lys_L_c    0.03   pchol_hs_c    3.3
# ser_L_c    0.017  2oxoadp_c     0.03
# phe_L_c    0.012  34hpp_c       0.011
# thr_L_c    0.012  urcan_c       0.01
# tyr_L_c    0.011  ac_c          0.005
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# trp_L_c    0.008
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

for i in range(369,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 780.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break
    
# =============================================================================
# 371 , ACACt2mB_MitoCore is a critical reaction
# testmodel2.summary()
# IN FLUXES          OUT FLUXES          OBJECTIVES
# -----------------  ------------------  --------------------
# fol_c      1e+03   co2_c        1e+03  OF_ATP_MitoCore  995
# h2o_c      1e+03   succ_c       1e+03
# so3_c    899       h_c        997
# ppa_c    477       biomass_c  995
# gln_L_c  292       so4_c      899
# o2_c      19.8     dhf_c      571
# hco3_c     1.71    thf_c      429
# hdca_c     1       4abut_c    240
# glc_D_c    0.9     orn_c      100
# lys_L_c    0.03    acac_c       4.66
# ser_L_c    0.017   pyr_c        1.86
# phe_L_c    0.012   2oxoadp_c    0.03
# thr_L_c    0.012   acald_c      0.012
# tyr_L_c    0.011   34hpp_c      0.011
# val_L_c    0.011
# asn_L_c    0.01
# his_L_c    0.01
# trp_L_c    0.008
# chol_c     0.0055
# met_L_c    0.005
# cys_L_c    0.001
# =============================================================================

testmodel2.reactions.get_by_id('ACACt2mB_MitoCore').bounds = -50.0, 1000.0
testmodel2.optimize()

for i in range(372,len(testmodel2.reactions)):       
    print(testmodel2.reactions[i].id, "bounds changed from", testmodel2.reactions[i].bounds, 
          "to", mymodel.reactions[i].bounds, ",", 'taken from', mymodel.reactions[i].id)
    testmodel2.reactions[i].bounds = mymodel.reactions[i].bounds 
    testmodel2.optimize()
    atp_flux = testmodel2.optimize().objective_value
    if atp_flux > 810.0:
        print(i, ',', testmodel2.reactions[i].id, 'is a critical reaction');break

testmodel = mymodel.copy()
testmodel.optimize()

for i in testmodel.reactions:
    if i.id != 'OF_ATP_MitoCore':
        if i.id != 'Biomass_MitoCore':
            i.lower_bound = 0.0

testmodel.optimize()

def change_critical_reaction_bounds_lb(maximum_atp_flux, model1 = testmodel, model2 = mymodel):
    #import time
    for i in range(1,len(model1.reactions)):       
        print(i, ',', model1.reactions[i].id, "bounds changed from", model1.reactions[i].bounds, 
              "to", model2.reactions[i].bounds, ",", 'taken from', model2.reactions[i].id)
        model1.reactions[i].bounds = model2.reactions[i].bounds 
        model1.optimize()
        #time.sleep(1)
        atp_flux = model1.optimize().objective_value
        if atp_flux > maximum_atp_flux:
            print(i, ',', model1.reactions[i].id, 'is a critical reaction')
            print(model1.reactions[i].id, "bounds changed from", model1.reactions[i].bounds, 
              "to -10.0, 1000.0")
            model1.reactions[i].bounds = -2.5, 1000.0
            atp_flux = model1.optimize().objective_value
            print('ATP flux = ', atp_flux)
            maximum_atp_flux = atp_flux + 0.5
        
list_of_critical_reactions_lb = [rec.id for rec in testmodel.reactions if rec.lower_bound == -2.5]        

testmodel = mymodel.copy()
testmodel.optimize()

for i in testmodel.reactions:
    if i.id != 'OF_ATP_MitoCore':
        if i.id != 'Biomass_MitoCore':
            i.upper_bound = 0.0

testmodel.optimize()

def change_critical_reaction_bounds_ub(maximum_atp_flux, model1 = testmodel, model2 = mymodel):
    import time
    for i in range(1,len(model1.reactions)):       
        print(i, ',', model1.reactions[i].id, "bounds changed from", model1.reactions[i].bounds, 
              "to", model2.reactions[i].bounds, ",", 'taken from', model2.reactions[i].id)
        model1.reactions[i].bounds = model2.reactions[i].bounds 
        model1.optimize()
        time.sleep(1)
        atp_flux = model1.optimize().objective_value
        if atp_flux > maximum_atp_flux:
            print(i, ',', model1.reactions[i].id, 'is a critical reaction')
            print(model1.reactions[i].id, "bounds changed from", model1.reactions[i].bounds, 
              'to', model1.reactions[i].lower_bound, ', 2.5')
            model1.reactions[i].upper_bound = 2.5
            atp_flux = model1.optimize().objective_value
            print('ATP flux = ', atp_flux)
            maximum_atp_flux = atp_flux + 0.5

testmodel.reactions.EX_o2_c.lower_bound = 0.0
testmodel.reactions.O2t.bounds = 0.0, 19.8
testmodel.reactions.O2tm.bounds = -1000.0, 1000.0

cobra.io.write_sbml_model(testmodel, "testmodel_29_05.xml")






