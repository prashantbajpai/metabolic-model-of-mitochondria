# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:16:45 2019

@author: prash
"""

#---------------------------------------------------------------------------------------------------
#load packages
#---------------------------------------------------------------------------------------------------
import cobra
import pandas as pd

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
mymodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/mitocore_updated_04_04.xml")
oldmodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/mitocore_updated_16_01.xml")
biggmodel = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/Recon3D_bigg.xml")
reconmodel = cobra.io.load_matlab_model("D:/My Ph.D/metabolic model of mitochondria/models/Recon3DModel_Dec2017.mat")


#---------------------------------------------------------------------------------------------------
#Of the metabolites needs to be tested, find which ones are now present in the model
#---------------------------------------------------------------------------------------------------
metabolites_to_test = pd.read_csv("D:/My Ph.D/metabolic model of mitochondria/metabolites_to_test.csv")
mymodel_all_metabolites = [met.id for met in mymodel.metabolites]
mets_present_in_mymodel = []
for i in metabolites_to_test.Recon_ID:    
    if i+"_c" in set(mymodel_all_metabolites):
        mets_present_in_mymodel.append(i+"_c")
    elif i+"_m" in set(mymodel_all_metabolites):
        mets_present_in_mymodel.append(i+"_m")

#---------------------------------------------------------------------------------------------------
#Find flux range of all metabolites
#---------------------------------------------------------------------------------------------------
        
from cobra.flux_analysis.variability import flux_variability_analysis


#too slow
#mymodel_bounds = [([rec.id for rec in mymodel.reactions], 
#                   [rec.lower_bound for rec in mymodel.reactions],
#                   [rec.upper_bound for rec in mymodel.reactions],
#                   [flux_variability_analysis(mymodel, reaction_list=rec.id, pfba_factor=1.1) for rec in mymodel.reactions])]
mymodel_flux_minimum = []
mymodel_flux_maximum = []
mymodel_upper_bound = []
mymodel_lower_bound = []
bar = myprogressbar(maxvalue=len(mymodel.reactions))
bar.start()
for i in range(len(mymodel.reactions)):
    mymodel_lower_bound.append(mymodel.reactions[i].lower_bound)
    mymodel_upper_bound.append(mymodel.reactions[i].upper_bound)
    bounds = flux_variability_analysis(mymodel, reaction_list=mymodel.reactions[i], pfba_factor=1.1)
    mymodel_flux_maximum.append(bounds.loc[mymodel.reactions[i].id, 'maximum'])
    mymodel_flux_minimum.append(bounds.loc[mymodel.reactions[i].id, 'minimum'])
    bar.update(i+1)
bar.finish()

mymodel_bounds = zip(mymodel_flux_minimum, mymodel_flux_maximum, mymodel_lower_bound,
                     mymodel_upper_bound)

import csv

with open("mitocore_bounds.csv", 'w', newline = '') as out:
    csv_out = csv.writer(out)
    csv_out.writerow(['Minimum Flux', 'Maximum Flux', 'Lower_bound', 'Upper_bound'])
    for rows in mymodel_bounds:
        csv_out.writerow(rows)

recon2_1 = cobra.io.read_sbml_model("D:/My Ph.D/metabolic model of mitochondria/models/MODEL1311110000_url.xml")

recon_flux_minimum = []
recon_flux_maximum = []
recon_upper_bound = []
recon_lower_bound = []
bar = myprogressbar(maxvalue=len(recon2_1.reactions))

bar.start()
for i in range(len(recon2_1.reactions)):
    recon_upper_bound.append(recon2_1.reactions[i].upper_bound)
    recon_lower_bound.append(recon2_1.reactions[i].lower_bound)
    bounds = flux_variability_analysis(recon2_1, reaction_list=recon2_1.reactions[i], pfba_factor=1.1)
    recon_flux_maximum.append(bounds.loc[recon2_1.reactions[i].id, 'maximum'])
    recon_flux_minimum.append(bounds.loc[recon2_1.reactions[i].id, 'minimum'])
    bar.update(i+1)
bar.finish()

recon2_bounds = zip(recon_flux_minimum, recon_flux_maximum, recon_lower_bound, recon_upper_bound)

with open("recon2_bounds.csv", 'w', newline = '') as out:
    csv_out = csv.writer(out)
    csv_out.writerow(['Minimum Flux', 'Maximum Flux', 'Lower_bound', 'Upper_bound'])
    for rows in recon2_bounds:
        csv_out.writerow(rows)

#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------        