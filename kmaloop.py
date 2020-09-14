#!/usr/bin/env python3
#Submits kma for isolates in specified directories. By default uses resfinder database
#190322 MAOS: removed "_L" requirement
import glob
import sys
import argparse
import os.path
import re
import os
from os import walk
import subprocess
from random import randint
import smtplib


parser = argparse.ArgumentParser(
	description='Looks for reads in subfolders to currentdir and submits them for KMA, default db is resfinder')
parser.add_argument("-p", "--partition", help="The queue that the job should be submitted to, be nice people!", default="daytime")
parser.add_argument("--db", help="The database that should be mapped against", default="/srv/data/DB/kma/resfinder_db")
parser.add_argument("-o", "--outputpattern", help="The overall structure of outputs", default=os.getcwd()+"/")
parser.add_argument("-d", "--directory", help="look in this directory", default=os.getcwd())
parser.add_argument("--rundir", help="Make and run in a new dir.", default=os.getcwd())
parser.add_argument("--filter", help="Filter qualifier for  which directories to go through, deprecated", default="")
parser.add_argument("--startfilter", help="Filter qualifier for how dirs start, deprecated", default="")
parser.add_argument("--filterlist", help="File containing the specific runs/folders that should be tested, deprecated")
parser.add_argument("--subdirs", help="Triggers behaviour to check in subdirectories.", action="store_true")
args = parser.parse_args()

notdonelist=[]
isolates={}
reflist=[]
testlist=[]
direction=["what?", "forward", "reverse"]
dir=os.getcwd()

def submitme (file, jobname):	#Function that checks links and submits the kma job
	try:
		info=re.search("([-\w]+)(_S\d+)(_L\d+)?(_R)(\d)(_001\.fastq\.gz)", file)
		if file.startswith("/"):
			path="/".join(file.split("/")[:-1])
			elements=file.split("/")
			if len(elements[-1]) > 2:
				file=elements[-1]
			else:
				file=elements[-2]
		if info and info.group(5)=="1":
			if not info.group(3):	#With or without lane number?
				os.system("".join(["sbatch -o /dev/null  --mem=1000 -c 2 -J ", jobname, " -p ", args.partition, " /srv/data/tools/git.repositories/SSI-scripts/maos/kmashell.sh '/tools/git.repositories/kma/kma -ipe ", path, "/", file, " ", path, "/", "".join([info.group(1), info.group(2), info.group(4), "2", info.group(6)]), " -t_db ", db, " -o ", args.outputpattern, info.group(1), "'"]))#-matrix 
			else:
				os.system("".join(["sbatch -o /dev/null  --mem=1000 -c 2 -J ", jobname, " -p ", args.partition, " /srv/data/tools/git.repositories/SSI-scripts/maos/kmashell.sh '/tools/git.repositories/kma/kma -ipe ", path, "/", file, " ", path, "/", "".join([info.group(1), info.group(2), info.group(3), info.group(4), "2", info.group(6)]), " -t_db ", db, " -o ", args.outputpattern, info.group(1), "'"]))#-matrix 
		else:
			info=re.search("(.+_R?)1(\.fastq\.gz)$", file)	#By no means a robust way of finding the forward read file, but should be fine in most cases.
			if info: 
				os.system("".join(["sbatch -o /dev/null --mem=1000  -J ", jobname, " -c 2 -p ", args.partition, " /srv/data/tools/git.repositories/SSI-scripts/maos/kmashell.sh '/tools/git.repositories/kma/kma -ipe ", path, "/", file, " ", path, "/", "".join([info.group(1), "2", info.group(2)]), " -t_db ", db, " -o ", args.outputpattern, info.group(1), "'"]))#-matrix 
				
	
	except:
		pass


#Insert test thingy, something somethin 
if not args.rundir:
	#Initiate some name for args.rundir. Based on args.directory.
	args.rundir=args.directory.split("/")[-1]
	if len(args.rundir) < 2:
		args.rundir=args.directory.split("/")[-2]

	
if args.rundir:
	#if not "/" in args.rundir:			#First we'll check that the new directory is in a specific enough place
	#	args.rundir="/".join(["/srv/data/FBI/tmp", args.rundir])
	#print(args.rundir)
	
	if not args.directory.startswith("/"):
		args.directory="/".join([os.getcwd(), args.directory])
	
	os.system("".join(["mkdir ", args.rundir]))	#Initiates the directory
	
	
	db=args.db
	if not os.path.isfile("".join([args.db, ".seq.b"])): 	#Check if the database is already formatted
		#Something something format the t_db file.
		db="/".join([args.rundir,"newdb"])
		os.system("".join(["/tools/git.repositories/kma/kma_index -i ", args.db, " -o ",db," >/dev/null 2>&1"]))
	os.chdir(args.rundir)						#Let's actually work in said directory	
	
	wgsrun=args.rundir.split("/")[-1]			#The last part of said dir should be the name of the run in question.
	args.outputpattern=os.getcwd()+"/"
	if len(wgsrun) < 2:							#Just to be safe regarding trailing / etc
		wgsrun=args.rundir.split("/")[-2]
	#os.system("".join(["ln -s ", args.directory, "/*gz ."])) #Make links to the reads ub said ru
	jobname="_".join(["kma", wgsrun])
	for file in os.listdir(args.directory):
		if "astq.gz" in file:
			#print(" ".join([file, jobname]))
			submitme("/".join([args.directory, file]), jobname)
	if args.db == "/srv/data/DB/kma/resfinder_db":	#If the database is resfinder, run ressummary
		
		os.system("".join(["sbatch -o /dev/null  --mem=1000 -c 2 -d singleton -J ", jobname, " -p ", args.partition, " /srv/data/tools/git.repositories/SSI-scripts/maos/ressummaryqcandkma.py -t -o ", args.rundir, "/ressummary.tsv"]))#-matrix 
	sys.exit()
	
if args.filterlist:
	filters=[]
	filterfile=open(args.filterlist, 'r')
	for filter in filterfile:
		filters.append(filter.rstrip())
	


if args.subdirs:

	for dir in os.listdir(args.directory):
		print(dir)
		if not dir.startswith(args.startfilter):
			continue
		if not args.filter in dir:
			continue
		if args.filterlist:
			filterstate="FilterMe"
			for filter in filters:
				if filter in dir:
					filterstate="LetMeIn"
					break
			if filterstate=="FilterMe":
				continue
		try:
			dir="/".join([args.directory, dir])
			jobname="Manualres"
			for file in os.listdir(dir):
				submitme(file, jobname)
				
		except:
			pass

		
		
#sbatch --mem=1000 -c 2 -p project /srv/data/tools/git.repositories/SSI-scripts/maos/kmashell.sh '/tools/git.repositories/kma/kma -ipe 2024-17-C1/2024-17-C1_R1.fastq.gz 2024-17-C1/C1_R2.fastq.gz -matrix -t_db /srv/data/DB/kma/resfinder_db -o /srv/data/FBI/ALZB/ecoli_subspecies/C1_R'