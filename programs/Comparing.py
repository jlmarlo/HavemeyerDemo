##creating input arguments for the scripts

import argparse
import json


parser = argparse.ArgumentParser(description='Comparing Annotation simularities between annotation programs')

parser.add_argument('annotation_file',type=str, help='Path to file containing a table of annotations from VEP, SNPEff, and ANNOVAR')
args = parser.parse_args()

##Files
table =open(args.annotation_file,'rt') #read in annotation table made by GATK
output = open('outputs/final/FinalReport.txt','wt') #output file
table.readline() #get rid of header
error = open('outputs/final/Errors.txt','wt') #error file
mismatch = open('outputs/final/MismatchLines.txt','wt') ##File to hold tables of mismatching variant classification and the counts of the number of times each disagreement happens
V_Amismatch = open('outputs/final/V_AMismatch.txt','wt') ##holds mismatches between Vep and Annovar MIGHT REMOVE OR NEED TO MAKE VEP SNPEFF VERSION
S_Amismatch = open('outputs/final/S_AMismatch.txt','wt') ##holds mismatches between SnpEff and Annovar SAME ABOVE
Somemismatch = open('outputs/final/3AnnotatorComparison.txt','wt') ##hold isntances of any disagreements
uneven = open('outputs/final/UniqueTranscripts.txt','wt') ##file to hold records of transcripts with no match
test = open('outputs/final/testoutput.txt','wt') ##file used for testing purposes
partial = open('outputs/final/partialMatches.txt','wt') ##file used to list partialmatches


linenumber=1 #set to one to account for header
##impact comparisons
same=0
diff=0
##Class comparisons between VEP and SNPEff
sameClass = 0
diffClass = 0
partClass = 0
##class comparisons between VEP and Annovar
V_Asame =0
V_Adiff =0
V_Apart = 0
##class comparisons between SnpEff and Annovar
S_Asame = 0
S_Adiff = 0
S_Apart = 0
##Class comparisons between all three annotators
AllSame=0
AllDiff = 0
AllPart = 0
SomeDiff = 0

simple =0
singleadditive = 0
simplelines =0
multitranscriptlines =0
multitranscript =0
additivemultitranscript =0

astericks = 0
## dictionary containing the common variant_Class and subclasses that will fall under the common class.
classdict = {
	"splicing_variant":['splicing','splicing_variant','splice_acceptor_variant','splice_donor_variant','splice_region_variant','splice_donor_5th_base_variant','splice_polypyrimidine_tract_variant'],
	"non_coding_transcript_variant":["non_coding_transcript_exon_variant","ncRNA_exonic",'ncRNA','non_coding_transcript_variant','non_coding_exon_variant','ncRNA_intronic','ncRNA_splicing','ncRNA_UTR5','ncRNA_UTR3'],
	"5_prime_UTR_variant":["5_prime_UTR_variant","UTR5",'5_prime_UTR_truncation','5_prime_UTR_premature_start_codon_gain_variant'],
	"3_prime_UTR_variant":["3_prime_UTR_variant","UTR3",'3_prime_UTR_truncation'],
	"intron_variant":['intronic_variant','intron','intronic','intron_variant','conserved_intron_variant'],
	"upstream_variant":['upstream_gene_variant','upstream'],
	"downstream_variant":["downstream_gene_variant","downstream"],
	"intergenic":['intergenic','intergenic_region','intergenic_variant','conserved_intergenic_variant'],
	"frameshift_variant":["frameshift_variant","frameshift_insertion","frameshift_deletion"],
	'stop_gained':['stop_gained','stopgain'],
	'stop_lost':['stop_lost','stoploss'],
	"inframe_insertion":['inframe_insertion','nonframeshift_insertion','disruptive_inframe_insertion','conservative_inframe_insertion'],
	"inframe_deletion":['inframe_deletion','disruptive_inframe_deletion','nonframeshift_deletion','conservative_inframe_deletion'],
	"nonsynonymous_variant":["missense_variant","nonsynonymous_SNV"],
	"synonymous_variant":["synonymous_variant","synonymous_SNV",'start_retained_variant','stop_retained_variant'],
	"start_lost":['start_lost'],
	"initiator_codon_variant":['initiator_codon_variant'],
	"coding_sequence_variant":["coding_sequence_variant"],
	"transcript_ablation":['transcript_ablation','feature_ablation'],
	"protein_altering_variant":['protein_altering_variant'],
	'exon_loss_variant':['exon_loss_variant'],
	'intragenic_variant':['intragenic_variant'],
	}

V_Amismatch.write('Line\tVEP\tANNOVAR')
S_Amismatch.write('Line\tSnpEff\tANNOVAR')
Somemismatch.write('Line\tSnpEff\tANNOVAR\tVEP')


EffimpactDict = {}
VEPimpactDict = {}
differences = {}
V_SClasses={}
V_AClasses={}
S_AClasses={}
partialDict={}
partialVADict={}
partialSADict = {}
bdict = {}
checkcount = 0


## Per variant comparison of annotations
for line in table: ##will loop through each line which represents a single variant


	VEPDict={} #initating dictionary for lines with multi transcripts
	EffDict={}
 
	infoList = line.split('\t') ##splitting columsn of info up.

	linenumber +=1 ##keeping track of what line number is being evaluated

	if '*' in infoList[3]: ##variants with * alleles are not evaluated by annotators so they will be removed
		error.write('Line number '+str(linenumber)+' removed due to * allele\n')
		astericks +=1 #keep track of how many lines are removed


	else:
		VEP = infoList[4] #get VEP annotation

		Eff = infoList[5] #get SNPEff annotation

		annovar = '|'.join(map(str,infoList[6:11])).strip() #grab and reconstruct annovar annotation


		## ANNOVAR annotations only have to be done one way per line
		ann_annotation = annovar.split('|') 	##split annotation into individual fields
		ann_area = ann_annotation[0] 		## the area of variant is stored in first field
		if '3b' in ann_area: 			##if '3b' is in field it is a "multitranscript" annotation
			multiarea= ann_area.split('\\x3b') 	##split on the 3b to get both parts of the annotation
			if 'exonic' in multiarea: 		##if exonic is one of the areas have to grab functional annotation
				multiarea.append(ann_annotation[3]) 	#add function to list of classes, would not include 2 areas and 1 function class
			commonlist = []
			
			##need to rename each item in the list of annotations to match common names
			for item in multiarea: #iterate through items in list
				#Change the variant classes into the common names
				for common,namePos in classdict.items():
					# for every common class in the dictionary, if the annotation is found in the list of children terms, change annotation
					# to parent term
					if item in namePos:
						commonName = common
						commonlist.append(commonName)
			ann_class = commonlist #set annotation class to the list of common class names
		
		##if there's not multitranscript, but area is exonic, variant_class is the function class. Change it to common name
		elif ann_area == 'exonic':
			ann_class=ann_annotation[3] #location of functional annotation
			for common,namePos in classdict.items():
				if ann_class in namePos:
					ann_class = common

		##if it's none of those things. The area is the variant class. Change to common name
		else:
			ann_class=ann_area
			for common,namePos in classdict.items():
				if ann_class in namePos:
					ann_class = common


		## For VEP and SnpEFf we evaluate annotation differently if there are multiple transcripts versus just one

		if ',' in VEP and ',' in Eff:  ##finding lines that have multiple transcripts

			TranSplitVEP = VEP.split(',') ##separating individual transcript predictions for VEP and SnpEff
			TranSplitEff = Eff.split(',')

			multitranscriptlines += 1
	
			##creating dictionaries for the multiple transcripts of VEP and SNPEff. Wills save transcript:class:impact
			for entryNum in range(len(TranSplitVEP)): #create dictionary for VEP transcripts
				CurrentVEP = TranSplitVEP[entryNum] ##grabs current transcript/annotation from list of transcript/annotations using index of list
				VEP_annotation = CurrentVEP.split('|') ##split annotation report into individual info fields
				VEPDict[VEP_annotation[6]]={'class':VEP_annotation[1],'impact':VEP_annotation[2]} # create dict entry to pair transcript with variant class and impact prediction
			for entryNum in range(len(TranSplitEff)): ##mirror of above but with snpeff
				CurrentEff = TranSplitEff[entryNum]
				Eff_annotation = CurrentEff.split('|')
				EffDict[Eff_annotation[6].split('.')[0]]={'class':Eff_annotation[1],'impact':Eff_annotation[2]} #transcript number contains version number that is removed for pairing


			##Determining if ther are transcripts not shared across VEP and SNPEff
			VEPunique = set(VEPDict.keys()) - set(EffDict.keys()) #Determining transcripts that one tool uses that the other tool does not
			Effunique = set(EffDict.keys()) - set(VEPDict.keys())

			##create reports of any nonshared transcripts
			if len(VEPunique) > 0: #creating records of transcripts that have no match
				for VU in VEPunique:
					print("VEP has a unique transcript:",VU,"with class", VEPDict[VU]['class'], "and impact",VEPDict[VU]['impact'],".\n",file=uneven)
					print(line,file=uneven)
			if len(Effunique) > 0: #creating records of transcripts that have no match
				for EU in Effunique:
					pass
					#print("SnpEff has a unique transcript:",EU,"with class", EffDict[EU]['class'], "and impact",EffDict[EU]['impact'],".\n",file=uneven)
			

			#### Making comparisons of annotations from transcripts that are found in both VEP and SnpEff
			SharedTranscripts = set(VEPDict.keys()) & set(EffDict.keys()) ##creating list of only transcripts present in both lists
			for transcript in SharedTranscripts:


				###Some annotations will have additive annotations. '&'
				###Annotations with & have to be split specially.
				##This top section will only be done on annotations with no &
				if '&' not in VEPDict[transcript]['class'] and '&' not in EffDict[transcript]['class']: ##passes on annotations that have & in either VEP or SNPEff

					
					#must normalize variant classification names between annotator programs
					for commonName, NamePos in classdict.items(): #iterating through key:value pairs in dictionary
						MatchFound = False
						if VEPDict[transcript]['class'] in NamePos: ##if the classification from VEP is in the list of possibilities for the current common name
							VEPDict[transcript]['common'] = commonName ##record common name in dictionary for this transcript
							MatchFound = True ##set match to true
							break ##stop cycle

					#this will check for any classes that don't have common names. This should not occur anymore
					if MatchFound == False:
						VEPDict[transcript]['common']=VEPDict[transcript]['class']
						print('No match found for VEP variant class',VEPDict[transcript]['class'],'Maintaing original classification',file=error)


					###Same as above but for SNPEff
					for commonName, NamePos in classdict.items():
						MatchFound = False
						if EffDict[transcript]['class'] in NamePos:
							EffDict[transcript]['common'] = commonName
							MatchFound = True
							break
					if MatchFound == False: 
						EffDict[transcript]['common']=EffDict[transcript]['class']
						print('No match found for SnpEff variant class', EffDict[transcript]['class'],'Maintaining original classification',file=error)
					
					##Following normalization record the impact and class corresponding to each tool of the current transcript
					VEP_impact = VEPDict[transcript]['impact']
					VEP_class = VEPDict[transcript]['common']
					Eff_impact = EffDict[transcript]['impact'] 
					Eff_class = EffDict[transcript]['common']

					##Counting up number of instances of impacts
					if VEP_impact in VEPimpactDict:
						VEPimpactDict[VEP_impact]+=1
					else:
						VEPimpactDict[VEP_impact]=1
				
					if Eff_impact in EffimpactDict:
						EffimpactDict[Eff_impact] +=1
					else:
						EffimpactDict[Eff_impact] =1	

					#Comparing Impact predictions and increase count of match or mismatched
					if VEP_impact == Eff_impact:
						same += 1 ##increases count of the amount of identical impact predictions
					else:
						diff += 1 ## increases count of the number of disagreeing impact predictions
						difference = VEP_impact+'_'+Eff_impact ##keeping track of which kind of mismatch is occuring LOW_MODERATE etc
						if difference in differences: ## if we've already seen one increase count of it
							differences[difference] += 1
						else:
							differences[difference] = 1 ##if not make new categories
					

					##comparing variant classes. Reminder these do not have '&'
					if VEP_class == Eff_class:
						sameClass += 1
					else:
						diffClass += 1
						hold = VEP_class+'_'+Eff_class
						if hold in V_SClasses:
							V_SClasses[hold] +=1
						else:
							V_SClasses[hold] = 1


					##Begins comparison of SnpEff and VEP with ANNOVAR
					##annotations are evaluated differently if ANNOVAR had multiple annotations
					if isinstance(ann_class,list): #if annovar annotation has multiple parts
						
						if VEP_class in ann_class: ##if the VEP classification is in the list of ANNOVAR annotations it's classed as a match
							V_Asame+=1
						else:
							V_Adiff +=1
						
						if Eff_class in ann_class:
							S_Asame +=1
						else:
							S_Adiff+=1
					
					else: #if annovar is a single prediction it's evaluated normally
						if VEP_class == ann_class:
							V_Asame +=1
						else:
							V_Adiff +=1
							hold = VEP_class+'_'+ann_class
							if hold in V_AClasses:
								V_AClasses[hold] +=1
							else:
								V_AClasses[hold] = 1

						if Eff_class == ann_class:
							S_Asame +=1
						else:
							S_Adiff +=1
							hold = Eff_class+'_'+ann_class
							if hold in S_AClasses:
								S_AClasses[hold] +=1
							else:
								S_AClasses[hold] = 1


					####evaluations of all three programs together. I need to work on this!
					if isinstance(ann_class,list):	
						if Eff_class in ann_class and VEP_class in ann_class and Eff_class == VEP_class:
							AllSame += 1
						elif Eff_class not in ann_class and VEP_class not in ann_class and Eff_class != VEP_class:
							AllDiff +=1
						else:
							SomeDiff +=1
					else:
						if Eff_class == ann_class and Eff_class == VEP_class and VEP_class == ann_class:
							AllSame +=1
						elif Eff_class != ann_class and Eff_class !=VEP_class and VEP_class != ann_class:
							AllDiff +=1
						else:
							SomeDiff +=1


		## ANNOTATIONS THAT ARE MULTITRANSCRIPT AND ADDITTIVE					
				else:
					VEP_impact = VEPDict[transcript]['impact']
					Eff_impact = EffDict[transcript]['impact']
				##initiating the percent match calcs. Might be unnecessary
					VSmatch =0
					VAmatch =0
					SAmatch = 0
				##
					##Counting number of impact classifications
					if VEP_impact in VEPimpactDict:
						VEPimpactDict[VEP_impact]+=1
					else:
						VEPimpactDict[VEP_impact]=1

					if Eff_impact in EffimpactDict:
						EffimpactDict[Eff_impact]+=1
					else:
						EffimpactDict[Eff_impact]=1
					##comparing effects
					if VEP_impact == Eff_impact:
						same += 1
					else:
						diff += 1
						difference = VEP_impact+'_'+Eff_impact
						if difference in differences:
							differences[difference] += 1
						else:
							differences[difference] = 1

					##Gather class information
					VEP_class = VEPDict[transcript]['class']
					Eff_class = EffDict[transcript]['class']
						
					##split the annotations on the & and create lists with all possibilities
					Vposs = VEP_class.split('&')
					Eposs = Eff_class.split('&')

					##Include loop to change all items in lists to common names
					commonVposs = []
					commonEposs = []
					for item in Vposs: #for all the annotations in VEP
						for common,namePos in classdict.items():
							if item in namePos: ##if the annotation is within a common category
								commonVposs.append(common) ##add common name to new list
					##add list to dictionary under 'common'
					VEPDict[transcript]['common'] = commonVposs

					for item in Eposs: ##same as above
						for common,namePos in classdict.items():
							if item in namePos:
								commonEposs.append(common)
					
					##add common name list to dictionary under 'common'
					EffDict[transcript]['common'] = commonEposs

					##find overlap of both lists of annotations. These are matches
					overlap = list(set(commonVposs) & set(commonEposs))
		
					##find values only found in one list of the other. These are mismatches
					Veponly = [x for x in commonVposs if x not in commonEposs] ##for every item in VEP if this item is not found in Eff add to Veponly
					Effonly = [x for x in commonEposs if x not in commonVposs]

					unique = len(Veponly)+len(Effonly) ##determining total number of classes unique to either list
					VSmatch = 100*(len(overlap)/(unique+len(overlap))) ##determine the percent of annotations that match

					if VSmatch == 100: ##if the match is 100 it will be considered a full match and added to same class count
						sameClass +=1
					elif VSmatch == 0: ##if it's 0% it's a mismatch
						diffClass += 1
					else: ##anything else is a partial match
						partClass += 1

						print(VSmatch,file=partial)
						print(overlap,Veponly,Effonly,sep='\t',file=partial)


						##keep track of what different levels of partial matches are being made
						if VSmatch in partialDict:
							partialDict[VSmatch] +=1
						else:
							partialDict[VSmatch] =1

					##evaluating ANNOVAR. Differently based on if the anovar annotation is a list or not. 
					if isinstance(ann_class,str):
						##if annovar annotation is a single item
						if ann_class in VEPDict[transcript]['common']:
							length = len(VEPDict[transcript]['common'])
							VAmatch = 100*(1/length)
							if VAmatch ==100:
								V_Asame += 1
							else:	
								V_Apart +=1
								if VAmatch in partialVADict:
									partialVADict[VAmatch] +=1
								else:
									partialVADict[VAmatch] =1
						else:
							
							V_Adiff +=1
						

						if ann_class in EffDict[transcript]['common']:
							length = len(EffDict[transcript]['common'])
							SAmatch = 100*(1/length)
							if SAmatch ==100:
								S_Asame +=1
							else:
								S_Apart +=1
								if SAmatch in partialSADict:
									partialSADict[SAmatch] += 1
								else:
									partialSADict[SAmatch] =1
						else:
							
							S_Adiff +=1
					else:
						VAoverlap = list(set(commonVposs) & set(ann_class))
						SAoverlap = list(set(commonEposs) & set(ann_class))
						
						VEPonly = [x for x in commonVposs if x not in ann_class]
						annonly = [x for x in ann_class if x not in commonVposs]
						unique = len(VEPonly)+len(annonly)
						VAmatch = 100*(len(VAoverlap)/(unique+len(VAoverlap)))

						if VAmatch == 100:
							V_Asame +=1
						elif VAmatch ==0:
							V_Adiff +=1
						else:
							V_Apart +=1
							if VAmatch in partialVADict:
								partialVADict[VAmatch] +=1
							else:
								partialVADict[VAmatch]=1

						Effonly = [x for x in commonEposs if x not in ann_class]
						annonly = [x for x in ann_class if x not in commonEposs]
						unique = len(Effonly)+len(annonly)
						SAmatch = 100*(len(SAoverlap)/(unique+len(SAoverlap)))
						if SAmatch ==100:
							S_Asame +=1
						elif SAmatch ==0:
							S_Adiff +=1
						else:
							S_Apart +=1
							if SAmatch in partialSADict:
								partialSADict[SAmatch] +=1
							else:
								partialSADict[SAmatch] =1
				###An evaluation of all three annotators. Need to evaluate this logic	
					if SAmatch == 100 and VAmatch ==100 and VSmatch == 100:
						AllSame +=1
					elif SAmatch ==0 and VAmatch ==0 and VSmatch ==0:
						AllDiff+=1
					else:
						SomeDiff +=1



		##COMPARISONS BEING MADE FOR VARIANTS WITH ONLY ONE TRANSCRIPT
		else:
	
			##Since there's only one transcript can jump straight into is
			VEP_annotation = VEP.split('|') #split VEP annotation into indv. info fields
			Eff_annotation = Eff.split('|') #split snpEff annotations into inv. info fields

			VEP_transcript= VEP_annotation[6]
			EFF_transcript= Eff_annotation[6].split('.')[0]

			VEP_class = VEP_annotation[1]
			VEP_impact = VEP_annotation[2]
			Eff_class = Eff_annotation[1]
			Eff_impact = Eff_annotation[2]
			
			#Counting up impacts for each tool
			if VEP_impact in VEPimpactDict:
				VEPimpactDict[VEP_impact] +=1
			else:
				VEPimpactDict[VEP_impact]=1
			
			if Eff_impact in EffimpactDict:
				EffimpactDict[Eff_impact] +=1
			else:
				EffimpactDict[Eff_impact]=1


			##There can still be additive annotations
			##Comparing differently if there is & or no &
			if '&' in VEP_class or '&' in Eff_class: ##If it is an additive annotation will need to split and perform partial matches

				###split on the ann and make a list of all pieces 
				Vposs = VEP_class.split('&')
				Eposs = Eff_class.split('&')
				commonVposs = []
				commonEposs = []

				### iterate and change names of all the pieces
				for item in Vposs:
					for common,namePos in classdict.items():
						if item in namePos:
							commonVposs.append(common)
				for item in Eposs:
					for common,namePos in classdict.items():
						if item in namePos:
							commonEposs.append(common)


				###perform partial match checks
				##matching vep and snpeff
				overlap = list(set(commonVposs) & set(commonEposs))
				Veponly = [x for x in commonVposs if x not in commonEposs]
				Effonly = [x for x in commonEposs if x not in commonVposs]
				unique = len(Veponly)+len(Effonly)
				VSmatch = 100*(len(overlap)/(unique+len(overlap)))
				if VSmatch == 100:
					sameClass += 1
				elif VSmatch == 0:
					diffClass += 1
				else:
					partClass +=1
					print(VSmatch,file=partial)
					print(overlap,Veponly,Effonly,sep='\t',file=partial)
					if VSmatch in partialDict:
						partialDict[VSmatch] +=1
					else:
						partialDict[VSmatch] = 1


				##matching vep and snpeff with annovar
				if isinstance(ann_class,str):
			
					if ann_class in commonVposs:
						VAmatch = 100*(1/len(commonVposs))
						if VAmatch ==100:
							V_Asame +=1
						else:
							V_Apart +=1
							if VAmatch in partialVADict:
								partialVADict[VAmatch] +=1
							else:
								partialVADict[VAmatch]=1
					else:
						V_Adiff +=1

					if ann_class in commonEposs:
						SAmatch = 100*(1/len(commonEposs))
						if SAmatch ==100:
							S_Asame +=1
						else:
							S_Apart +=1
							if SAmatch in partialSADict:
								partialSADict[SAmatch]+=1
							else:
								partialSADict[SAmatch]=1
					else:
						S_Adiff+=1
				else:
					VAoverlap = list(set(commonVposs) & set(ann_class))
					SAoverlap = list(set(commonEposs) & set(ann_class))
					print(set(ann_class),set(commonEposs))
					VEPonly = [x for x in commonVposs if x not in ann_class]
					annonly = [x for x in ann_class if x not in commonVposs]
					unique = len(VEPonly)+len(annonly)
					VAmatch = 100*(len(VAoverlap)/(unique+len(VAoverlap)))

					if VAmatch ==100:
						V_Asame +=1
					elif VAmatch ==0:
						V_Adiff +=1
					else:
						V_Apart +=1
						if VAmatch in partialVADict:
							partialVADict[VAmatch] +=1
						else:
							partialVADict[VAmatch] =1


					Effonly = [x for x in commonEposs if x not in ann_class]
					annonly = [x for x in ann_class if x not in commonEposs]
					unique = len(Effonly)+len(annonly)
					SAmatch = 100*(len(SAoverlap)/(unique+len(SAoverlap)))
				if SAmatch ==100:
					S_Asame +=1
				elif SAmatch ==0:
					S_Adiff += 1
				else:
					S_Apart +=1
					if SAmatch in partialSADict:
						partialSADict[SAmatch]+=1
					else:
						partialSADict[SAmatch] =1
		


				##trying to evaluate all three programs together. will need to be checked	
				if SAmatch ==100 and VAmatch==100 and VSmatch==100:
					AllSame +=1
				elif SAmatch ==0 and VAmatch ==0 and VSmatch ==0:
					AllDiff +=1
				else:
					SomeDiff +=1	
					


			else: ##if there aren't any & can proceed as normal
				##Changing the names of all the variant classes for all programs
				for commonName, NamePos in classdict.items():
					if VEP_annotation[1] in NamePos:
						VEP_class = commonName
					if Eff_annotation[1] in NamePos:
						Eff_class = commonName


					
				##Comparing VEP and Eff class
				if VEP_class == Eff_class:
					sameClass +=1
				else:
					diffClass +=1
					hold = VEP_class+'_'+Eff_class
					if hold in V_SClasses:
						V_SClasses[hold] +=1
					else:
						V_SClasses[hold]=1
					#print(str(linenumber),'NA',VEP_class,Eff_class,sep='\t',file=mismatch)
			
				##Comparing VEP and Ann
				if VEP_class == ann_class:
					V_Asame +=1
				else:
					V_Adiff +=1
					print('VEP',VEP_class,ann_class,file=test)
					hold = VEP_class+'_'+ann_class
					if hold in V_AClasses:
						V_AClasses[hold]+=1
					else:
						V_AClasses[hold]=1
					print(linenumber,VEP_class,ann_class,sep='\t',file=V_Amismatch)

				##Comparing Eff and Ann
				if Eff_class == ann_class:
					S_Asame +=1
				else:
					S_Adiff +=1
					print('Eff',Eff_class,ann_class,file=test)
					hold = Eff_class+'_'+ann_class
					if hold in S_AClasses:
						S_AClasses[hold] +=1
					else:
						S_AClasses[hold]=1
					print(linenumber,Eff_class,ann_class,sep='\t',file=S_Amismatch)


				##Comparing all three
				if Eff_class == ann_class == VEP_class:
					AllSame +=1
				else:
					SomeDiff +=1
					print(linenumber,Eff_class,ann_class,VEP_class,sep='\t',file=Somemismatch)
			##impact
				
			if VEP_impact == Eff_impact:
				same +=1
			else: 
				diff +=1
				#print(str(linenumber),VEP_transcript,VEP_impact,Eff_impact,sep='\t',file=mismatch)
				difference = VEP_impact+'_'+Eff_impact
				if difference in differences:
					differences[difference] +=1 
				else:
					differences[difference]=1





######Creating report file#######
ImpactMatchPer = 100*(same/(same+diff))
ImpactMisMatchPer = 100*(diff/(same+diff))


ClassMatchPer = 100*(sameClass/(sameClass+diffClass+partClass))
ClassMisMatchPer = 100*(diffClass/(sameClass+diffClass+partClass))
ClassPartMatch = 100*(partClass/(sameClass+diffClass+partClass))

ClassMatchV_A = 100*(V_Asame/(V_Asame+V_Adiff+V_Apart))
ClassPartMatchV_A = 100*(V_Apart/(V_Asame+V_Adiff+V_Apart))
ClassMismatchV_A = 100*(V_Adiff/(V_Asame+V_Adiff+V_Apart))

ClassMatchS_A = 100*(S_Asame/(S_Asame+S_Adiff+S_Apart))
ClassPartMatchS_A = 100*(S_Apart/(S_Asame+S_Adiff+S_Apart))
ClassMismatchS_A = 100*(S_Adiff/(S_Asame+S_Adiff+S_Apart))

ClassMatchAll = 100*(AllSame/(AllSame+SomeDiff+AllDiff))
ClassPartMatchAll = 100*(SomeDiff/(AllSame+SomeDiff+AllDiff))
ClassMismatchAll = 100*(AllDiff/(AllSame+SomeDiff+AllDiff))

print('There were',astericks,'number of lines removed for invalid alleles',file=output)

output.write(str(AllSame)+' '+str(AllDiff)+' '+str(SomeDiff))
output.write('There were ' +str(same)+ ' matching impact predictions ('+ str(ImpactMatchPer)+ '%).\nThere were ' + str(diff) + ' different impact predictions ('+str(ImpactMisMatchPer)+'%).\n')
output.write('VEP_SnpEff\nThere were ' +str(sameClass)+ ' matching variant class predictions ('+str(ClassMatchPer)+ '%).\nThere were '+str(partClass)+ ' partial variant class predictions ('+str(ClassPartMatch)+'\n'+'There were ' + str(diffClass) + ' different variant class predictions ('+ str(ClassMisMatchPer)+ '%).\n\n')


output.write('VEP_ANNOVAR\nThere were '+str(V_Asame)+ ' matching variant class predictions ('+str(ClassMatchV_A)+'%).\n'+'There were '+str(V_Apart)+' partial variant class predictions ('+str(ClassPartMatchV_A)+'%).\n'+'There were ' + str(V_Adiff) + 'difference variant class predictions ('+str(ClassMismatchV_A)+'%).\n\n')


output.write('SnpEff_ANNOVAR\nThere were '+str(S_Asame)+ ' matching variant class predictions ('+str(ClassMatchS_A)+'%).\n'+'There were'+str(S_Apart)+' partial variant class predictions ('+str(ClassPartMatchS_A)+'%).\n'+'There were ' + str(S_Adiff) + 'different variant class predictions ('+str(ClassMismatchS_A)+'%).\n\n')


output.write('VEP_ANNOVAR_SnpEff\nThere were '+str(AllSame)+ ' matching variant class predictions ('+str(ClassMatchAll)+'%).\nThere were ' + str(SomeDiff) + 'instances of partial matching ('+str(ClassPartMatchAll)+'%).\n There were '+ str(AllDiff)+' completely different  variant class predictions ('+str(ClassMismatchAll)+'%).\n\n')

print('VEP Impact\tSnpEff Impact\tNumber of Occurences',file=output)
for key in differences.keys():
	first,second = key.split('_')
	occurence = differences[key]
	occurencepercent = 100*(occurence/sum(differences.values()))
	print(first,second,occurence,occurencepercent,sep='\t',file=output)
print('Total Number of Impact Differences:',sum(differences.values()),file=output)

print('VEP Impact\tCount',file=output)
for key,value in VEPimpactDict.items():
	print(key,':',value,file=output)
print('Total Impact predictions by VEP:',sum(VEPimpactDict.values()),file=output)


print('Eff Impact\tCount',file=output)
for key,value in EffimpactDict.items():
	print(key,':',value,file=output)
print('Total Impact predictions by SnpEff:',sum(EffimpactDict.values()),file=output)

#for key in Classes.keys():
#	print(key,Classes[key],file=output)
print("\n\nPartial States VEP_SNPEff",file=output)

sortedPartialVS = dict(sorted(partialDict.items()))
for key,value in sortedPartialVS.items():
	percentValue = 100*(value/sum(sortedPartialVS.values()))
	print(key,":",value,percentValue,file=output)
print('Total Number of Partial Matches:',sum(sortedPartialVS.values()),file=output)




print("\n\nPartial Stats VEP_Annovar",file=output)
sortedPartialVA = dict(sorted(partialVADict.items()))
for key,value in sortedPartialVA.items():
	percentValue = 100*(value/sum(sortedPartialVA.values()))
	print(key,":",value,percentValue,file=output)
print('Total Number of Partial Matches:',sum(sortedPartialVA.values()),file=output)



print("\n\nPartial Stats SNPEff_Annovar",file=output)
sortedPartialSA = dict(sorted(partialSADict.items()))
for key,value in sortedPartialSA.items():
	percentValue = 100*(value/sum(sortedPartialSA.values())) 
	print(key,':',value,percentValue,file=output)
print('Total Number of Partial Matches:',sum(sortedPartialSA.values()),file=output)


print('VEP','SNPEff','Count',sep='\t',file=mismatch)
for key,value in V_SClasses.items():
	percentValue = 100*(value/sum(V_SClasses.values()))
	print(key,value,percentValue,sep='\t',file=mismatch)
print("Total Number of Mismatches:",sum(V_SClasses.values()),file=mismatch)

print('\n\n\n\n',file=mismatch)
print('VEP','ANNOVAR','Count',sep='\t',file=mismatch)
for key,value in V_AClasses.items():
	percentValue = 100*(value/sum(V_AClasses.values()))
	print(key,value,percentValue,sep='\t',file=mismatch)
print("Total Number of Mismatches:",sum(V_AClasses.values()),file=mismatch)
print('\n\n\n\n',file=mismatch)
print('SNPEff','ANNOVAR','Count',sep='\t',file=mismatch)
for key,value in S_AClasses.items():
	percentValue = 100*(value/sum(S_AClasses.values()))
	print(key,value,percentValue,sep='\t',file=mismatch)
print('Total Number of Mismatches:',sum(S_AClasses.values()),file=mismatch)


#print('Total Difference from VEP:',VEP_part_diff,"Partial Match with VEP:",VEP_part_same,"Total Difference from Eff:",Eff_part_diff,"Partial Match with Eff:",Eff_part_same,"Matched Nothing:",totalWhiff,sep='\n',file=test)
#print('Total lines with multiple transcripts',multitranscriptlines,'Comparisons from multiple transcript lines with annotations that are simple',multitranscript,'Multitranscript additive comparisons',additivemultitranscript,'Total lines with single transcripts',simplelines,'Single transcript non additive comparisons',simple,'Single transcript with additive annotations',singleadditive,sep='\n',file=output)

table.close()
output.close()
error.close()
mismatch.close()
uneven.close()
test.close()
V_Amismatch.close()
S_Amismatch.close()
Somemismatch.close()
partial.close()
