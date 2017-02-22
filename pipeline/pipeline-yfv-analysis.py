#################################################################
#################################################################
###############  ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys
import pandas as pd
import numpy as np
import rpy2.robjects as robjects
import pandas.rpy.common as com

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support as S
import PipelineYfvAnalysis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
seriesMatrixFile = 'rawdata.dir/expression/GSE51972_series_matrix.txt'
clinicalAnnotationFile = 'rawdata.dir/annotations/YFV_clinicalparameters.xls'
platformAnnotationFile = 'rawdata.dir/annotations/GPL3535-10024.txt'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-yfv-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Annotation Data
#######################################################
#######################################################

#############################################
########## 1. Sample Annotation Data
#############################################

@follows(mkdir('f1-annotation_data.dir'))

@merge([seriesMatrixFile,
		clinicalAnnotationFile],
	   'f1-annotation_data.dir/yfv-sample_annotations.txt')

def processSampleAnnotations(infiles, outfile):

	# Open file
	with open(seriesMatrixFile, 'r') as openfile:

	    # Read data in list
	    dataList = [x.replace('!', '').replace('"', '').strip() for x in openfile.readlines() if x[0] == '!']

	# Convert to dict
	dataDict = {x.split('\t')[0]:x.split('\t')[1:] for x in dataList}

	# Select keys
	selectedKeys = ['Sample_characteristics_ch1',
	                'Sample_geo_accession',
	                'Sample_source_name_ch1',
	                'Sample_title']

	# Get subset
	dataDictSubset = {x:dataDict[x] for x in dataDict.keys() if x in selectedKeys}

	# Convert to dataframe
	annotationDataframe = pd.DataFrame(dataDictSubset)

	# Fix columns
	annotationDataframe['age'] = [x.replace('age: ', '') for x in annotationDataframe['Sample_characteristics_ch1']]
	annotationDataframe['Sample_title'] = [x.replace(' ', '-') for x in annotationDataframe['Sample_title']]

	# Add columns
	annotationDataframe['treatment'] = [x.split('_')[1] for x in annotationDataframe['Sample_title']]
	annotationDataframe['timepoint'] = [x.split('_')[2] for x in annotationDataframe['Sample_title']]
	annotationDataframe['replicate'] = [x.split('_')[3] for x in annotationDataframe['Sample_title']]
	annotationDataframe['DPI'] = [int(x.replace('day', '')) for x in annotationDataframe['timepoint']]

	# Read blood values dataframe
	bloodMeasurementDataframe = pd.read_excel(clinicalAnnotationFile, 0)

	# Read sample matching dataframe
	sampleMatchingDataframe = pd.read_excel(clinicalAnnotationFile, 1)

	# Add column
	sampleMatchingDataframe['replicate'] = ['rep%(x)s' % locals() for x in sampleMatchingDataframe['Replicate# for microarray study']]

	# Merge
	mergedBloodMeasurementDataframe = bloodMeasurementDataframe.merge(sampleMatchingDataframe, on='Animal ID')

	# Merge results
	mergedDataframe = annotationDataframe.merge(mergedBloodMeasurementDataframe, on=['replicate', 'DPI'])

	# Rename columns
	mergedDataframe.rename(columns={'Day of Euthanasia':'day_of_euthanasia', 'Animal ID':'animal_ID', 'Viral Loads':'viral_loads'}, inplace=True)

	# Get specified columns
	mergedDataframe.drop(['Sample_characteristics_ch1', 'Replicate# for microarray study'], axis=1, inplace=True)

	# Write file
	mergedDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Gene Annotation Data
#############################################

@files(platformAnnotationFile,
	   'f1-annotation_data.dir/GPL3535-gene_annotations.txt')

def getGeneAnnotations(infile, outfile):

	# Read dataframe
	geneAnnotationDataframe = pd.read_table(infile, comment='#')

	# Get columns
	columnDict = {'ID': 'probe_id', 'Gene Symbol': 'gene_symbol', 'ENTREZ_GENE_ID': 'entrez_gene_id', 'Species Scientific Name': 'species'}

	# Get subset
	geneAnnotationDataframe = geneAnnotationDataframe[columnDict.keys()].rename(columns=columnDict)

	# Save
	geneAnnotationDataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S2. Expression Data
#######################################################
#######################################################

#############################################
########## 3. Gene Expression Data
#############################################

@follows(mkdir('f2-expression_data.dir'))

@files([seriesMatrixFile,
		getGeneAnnotations],
	   'f2-expression_data.dir/yfv-gene_expression_matrix.txt')

def getGeneExpressionData(infiles, outfile):

	# Split infiles
	seriesMatrixFile, geneAnnotationFile = infiles

	# Get expression dataframe
	expressionDataframe = pd.read_table(seriesMatrixFile, comment='!', index_col='ID_REF')

	# Get annotation dataframe
	geneAnnotationDataframe = pd.read_table(geneAnnotationFile, index_col='probe_id')

	# Get probe to gene conversion dict
	probe2gene = geneAnnotationDataframe.loc[geneAnnotationDataframe['species'] == 'Rhesus macaque', ['gene_symbol']].dropna().to_dict()['gene_symbol']

	# Get gene symbols
	expressionDataframe['gene_symbol'] = [probe2gene[x] if x in probe2gene.keys() else np.nan for x in expressionDataframe.index]

	# Remove NaN
	expressionDataframeFiltered = expressionDataframe.dropna()

	# Group by mean
	expressionDataframeGrouped = expressionDataframeFiltered.groupby('gene_symbol').mean()

	# Save
	expressionDataframeGrouped.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################






#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
