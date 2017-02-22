#################################################################
#################################################################
############### Yellow Fever Analysis ###########################
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
humanHomologyFile = 'rawdata.dir/annotations/macaque_human_homologs.txt'

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

@follows(mkdir('f1-annotation.dir'))

@merge([seriesMatrixFile,
		clinicalAnnotationFile],
	   'f1-annotation.dir/yfv-sample_annotations.txt')

def getSampleAnnotations(infiles, outfile):

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
	   'f1-annotation.dir/GPL3535-gene_annotations.txt')

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

@follows(mkdir('f2-raw_expression.dir'))

@files([seriesMatrixFile,
		getGeneAnnotations],
	   'f2-raw_expression.dir/yfv-gene_expression.txt')

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
########## S3. Run Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Run Characteristic Direction
#############################################

@follows(mkdir('f3-differential_expression.dir/cd'))


@transform(getGeneExpressionData,
		   regex(r'.*/(.*)_expression.txt'),
		   add_inputs(getSampleAnnotations),
		   r'f3-differential_expression.dir/cd/\1_differential_expression.txt')

def runCharacteristicDirection(infiles, outfile):

	# Split infiles
	expressionFile, annotationFile = infiles

	# Read expression data
	expressionDataframe = pd.read_table(expressionFile, index_col='gene_symbol')

	# Read annotation data
	annotationDataframe = pd.read_table(annotationFile, index_col='Sample_geo_accession')

	# Reindex annotation dataframe
	annotationDataframeReindex = annotationDataframe.reset_index().set_index(['treatment','timepoint'])

	# Get unique treatments and timepoints
	treatments = set(annotationDataframe['treatment'])
	timepoints = set(annotationDataframe['timepoint'])

	# Get sample dictionary
	sampleDict = {'_'.join([treatment, timepoint]): annotationDataframeReindex.loc[(treatment, timepoint), 'Sample_geo_accession'].tolist() for timepoint in timepoints for treatment in treatments}

	# Initialize empty dataframe
	resultDataframe = pd.DataFrame()

	# Loop through comparisons
	for comparisons in [['17D-vaccinated_day0', '17D-vaccinated_day3'], ['YFVwt-infected_day0', 'YFVwt-infected_day3'], ['17D-vaccinated_day0', 'YFVwt-infected_day0'], ['17D-vaccinated_day3', 'YFVwt-infected_day3']]:
	    
		# Get columns
		controlColumns, experimentColumns = [sampleDict[x] for x in comparisons]

		# Run characteristic direction
		cdResults = r.runCharacteristicDirection(com.convert_to_r_dataframe(expressionDataframe), experimentColumns, controlColumns, 0.01)

		# Convert to dataframe
		cdDataframe = com.convert_robj(cdResults).reset_index()

		# Add comparison column
		cdDataframe['comparison'] = 'v'.join(comparisons)

		# Append
		resultDataframe = pd.concat([resultDataframe, cdDataframe])

	# Pivot
	resultDataframeCast = resultDataframe.pivot(index='index', columns='comparison', values='CD')

	# Save
	resultDataframeCast.to_csv(outfile, sep='\t', index_label='gene_symbol')


#######################################################
#######################################################
########## S4. Enrichment Analysis
#######################################################
#######################################################

#############################################
########## 1. Submit Genesets 
#############################################

@follows(mkdir('f4-enrichr.dir'))

@transform(runCharacteristicDirection,
		   regex(r'.*/(.*)differential_expression.txt'),
		   r'f4-enrichr.dir/\1enrichr_links.txt')

def submitEnrichrGenesets(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile, index_col='gene_symbol').fillna(0)

	# Rewrite index
	cdDataframe.index = [x.split(' /// ')[0] for x in cdDataframe.index]

	# Initialize link dataframe
	resultDataframe = pd.DataFrame()

	# Loop through timepoints
	for comparison in cdDataframe.columns:

	    # Get Enrichr links
	    enrichrLinkDataframe = S.uploadToEnrichr(cdDataframe, comparison)

	    # Add timepoint label
	    enrichrLinkDataframe['comparison'] = comparison

	    # Concatenate
	    resultDataframe = pd.concat([resultDataframe, enrichrLinkDataframe])

	# Save data
	resultDataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Get Enrichment results
#############################################

@transform(submitEnrichrGenesets,
		   regex(r'.*/(.*)links.txt'),
		   r'f4-enrichr.dir/\1results.txt')

def getEnrichrResults(infile, outfile):

	# Read infile
	enrichrLinkDataframe = pd.read_table(infile, index_col=['geneset','comparison'])

	# Initialize result dataframe
	resultDataframe = pd.DataFrame()

	# Set libraries
	libraries = ['ChEA_2016', 'KEGG_2016', 'GO_Biological_Process_2015', 'GO_Cellular_Component_2015', 'GO_Molecular_Function_2015', 'VirusMINT']

	# Loop through timepoints, genesets and libraries
	for geneset in enrichrLinkDataframe.index.levels[0]:
	    for comparison in enrichrLinkDataframe.index.levels[1]:
	        for library in libraries:

	            # Get enrichment results
	            enrichmentResultDataframe = S.getEnrichmentResults(enrichrLinkDataframe.loc[(geneset, comparison), 'userListId'], library)

	            # Add labels
	            enrichmentResultDataframe['comparison'] = comparison
	            enrichmentResultDataframe['geneset'] = geneset
	            enrichmentResultDataframe['library'] = library

	            # Concatenate
	            resultDataframe = pd.concat([resultDataframe, enrichmentResultDataframe])

    # Write file
	resultDataframe.to_csv(outfile, sep='\t', index=False)


#######################################################
#######################################################
########## S9. L1000CDS2 Analysis
#######################################################
#######################################################

#############################################
########## 1. Submit analysis
#############################################

@follows(mkdir('f5-l1000cds2.dir'))

@transform(runCharacteristicDirection,
		   regex(r'.*/(.*)differential_expression.txt'),
		   r'f5-l1000cds2.dir/\1l1000cds2_links.txt')

def runL1000CDS2(infile, outfile):

	# Read infile
	cdDataframe = pd.read_table(infile, index_col='gene_symbol').fillna(0)

	# Rewrite index
	cdDataframe.index = [x.split(' /// ')[0] for x in cdDataframe.index]

	# Initialize dataframes
	linkDataframe = pd.DataFrame()
	signatureDataframe = pd.DataFrame()

	# Loop through timepoints
	for comparison in cdDataframe.columns:

	    # Run L1000CDS2
	    resultDict = S.getL1000CDS2Results(cdDataframe, comparison)

	    # Add comparison labels
	    resultDict['links']['comparison'] = comparison
	    resultDict['signatures']['comparison'] = comparison

	    # Append dataframes
	    linkDataframe = pd.concat([linkDataframe, resultDict['links']])
	    signatureDataframe = pd.concat([signatureDataframe, resultDict['signatures']])

	# Write files
	linkDataframe.to_csv(outfile, sep='\t', index=False)
	signatureDataframe.to_csv(outfile.replace('links', 'signatures'), sep='\t', index=False)


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
