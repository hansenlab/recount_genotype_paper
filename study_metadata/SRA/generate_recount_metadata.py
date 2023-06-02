#!/usr/bin/env python3


import os, sys, warnings, pathlib
import pandas as pd


def searchDirectoryGenerator(searchDir, level=1, fixed=False):
	"""Returns a generator of 3-tuple (root, dirs, files) no deeper than a specified directory level.
	Inspired from https://stackoverflow.com/questions/229186/os-walk-without-digging-into-directories-below
	Keyword arguments:
	searchDir -- the base directory to start searching
	level -- ceiling number of directory levels go down searching
	fixed -- boolean of whether we yield only directory at only level-depth. 
			 default is to return at most level-depth.
    """
	searchDir = searchDir.rstrip(os.path.sep)
	assert os.path.isdir(searchDir), "Input is not a directory."
	searchDirLevel = searchDir.count(os.path.sep) #count number of directories in searchDir via "/"
	for root, dirs, files in os.walk(searchDir):
		if not fixed:
			yield root, dirs, files
		currentDirLevel = root.count(os.path.sep)
		if currentDirLevel > searchDirLevel + level:
			del dirs[:]
		if fixed and searchDirLevel + level == currentDirLevel:
			yield root, dirs, files

def findFile(searchDir, sampleName, projectName, endsWith):
	"""Find a file from `searchDir` that ends with `endsWith` and contains the `sampleName` and `projectName`"""
	for root, dirs, files in os.walk(searchDir):
		for f in files:
			if f.endswith(endsWith):
				if f.find(sampleName) == -1:
					warnings.warn("Found file", f, "that doesn't belong within sample name: ", sampleName)
				if f.find(projectName) == -1:
					warnings.warn("Found file", f, "that doesn't belong within project name: ", projectName)
				return os.path.join(root, f)
	return None

def generateProjectMetadata(project_ID, SEARCH_DIR):
	project_ID = project_ID.rstrip().upper()
	assert project_ID.startswith("SRP") or project_ID.startswith("ERP") or \
		project_ID.startswith("DRP"), "Project ID must start with SRP or ERP or DRP."
	
	metadata = {}
	#Search for the project.
	print("Searching for:", project_ID)
	foundProjectDir = ""
	project_generator = searchDirectoryGenerator(SEARCH_DIR, level=2)
	for root, dirs, files in project_generator:
		for name in dirs:
			if name.strip().upper() == project_ID:
				foundProjectDir = os.path.join(root, name)
				print(project_ID, "found! At:", foundProjectDir)
				break
	if foundProjectDir == "":
		print("Could not find", project_ID, "under", SEARCH_DIR)
		return None

	#Search for samples.
	sample_generator = searchDirectoryGenerator(foundProjectDir, level=1, fixed=True)
	for root, dirs, files in sample_generator:
		for d in dirs:
			sampleName = d.strip().upper()
			if (project_ID.startswith("SRP") and sampleName.startswith("SRR")) or \
			   (project_ID.startswith("ERP") and sampleName.startswith("ERR")) or \
			   (project_ID.startswith("DRP") and sampleName.startswith("DRR")):
				metadata[sampleName] = {'total': findFile(os.path.join(root, d), sampleName, project_ID, ".all.bw"), 
										'alt': findFile(os.path.join(root, d), sampleName, project_ID, ".bamcount_nonref.csv.zst")}
			else:
				warnings.warn("Warning: Found sample folder that's not SRR or ERR or DRR:", os.path.join(root, d))

	#Reformat, save, exit.
	df = pd.DataFrame.from_dict(metadata, orient='index')
	df.index.name = 'sample_id'
	df.reset_index(inplace=True)
	df['study'] = project_ID
	return df



if __name__ == "__main__":
	SEARCH_DIR = "/dcl02/lieber/ajaffe/recount-pump/human_tranche_backups" # where SRA data is stored.
	OUTPUT_FOLDER = "/dcs04/hansen/data/recount_genotype/pipeline/metadata/"

	recount_projects = pd.read_csv("recount3_noGTEx_noTCGA.csv") #all human SRA from recount (excluding TCGA and GTEx)

	progress = 0
	for project in recount_projects['project']:
		result = generateProjectMetadata(project, SEARCH_DIR)
		result.to_csv(OUTPUT_FOLDER + project + ".csv", index=False)

		progress += 1
		print(progress)

	print("Finished!")
	quit()
