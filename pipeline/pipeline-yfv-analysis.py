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

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-yfv-analysis.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Annotation Data
#############################################

@follows(mkdir('f1-processed_data.dir'))

@files('rawdata.dir/GSE51972_series_matrix.txt',
	   'f1-processed_data.dir/yfv-sample_annotations.txt')

def processSampleAnnotations(infile, outfile):

	# Open file
	with open(infile, 'r') as openfile:
	    
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
	annotationDataframe.drop('Sample_characteristics_ch1', axis=1, inplace=True)

	# Add columns
	annotationDataframe['treatment'] = [x.split('_')[1] for x in annotationDataframe['Sample_title']]
	annotationDataframe['timepoint'] = [x.split('_')[2] for x in annotationDataframe['Sample_title']]
	annotationDataframe['replicate'] = [x.split('_')[3] for x in annotationDataframe['Sample_title']]

	# Write file
	annotationDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Expression Data
#############################################

@files('rawdata.dir/GSE51972_series_matrix.txt',
	   'f1-processed_data.dir/yfv-expression_matrix.txt')

def processExpressionData(infile, outfile):

	# Read data
	seriesDataframe = pd.read_table(infile, comment='!')

	# Save data
	seriesDataframe.to_csv(outfile, sep='\t', index=False)

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
