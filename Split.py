#Name:          Mark TK Cheng
#Email:         NA
#Date:          06th Sept 2021
#Function(s):   Splits multiallelels output of snp-sites
#Inspiration:   https://www.biostars.org/p/302940/

import sys
import os

if not len(sys.argv) ==3:
    print ("\nError:\tincorrect number of command-line arguments")
    print ("Syntax:\tSplit.py [Input VCF] [Output VCF]\n")
    sys.exit()

if sys.argv[1] == sys.argv[2]:
    #boolOverwrite = raw_input("\nInput file is the same as the output file - overwrite input file (sim/nao)?\n")
    print ("Error:\tInput file is the same as the output file - choose a different output file\n")
    sys.exit()

#File input 
fileInput = open(sys.argv[1],"r")
fileOutput = open(sys.argv[2],"w")


#Loop through each line in the input file, and split multiallelic sites
print ("Splitting multi-allelic sites...")
for strLine in fileInput.readlines():
    #Strip the endline character from each input line
    print(strLine)
    strLine = strLine.strip("\"")
    #The '#' character in VCF format indicates that the line is a header. Ignore thse and just output to the new file
    if strLine.startswith("#"):
        fileOutput.write(strLine)
    else:
        #split the tab-delimited line into an array
        strArray = [splits for splits in strLine.split("\t") if splits != ""]
        strArray[-1] = strArray[-1].strip() #strips /n in the last string in the array
        
        #check first it it's multiallelic
        if "," in strArray[4]:
            strVars = [splits for splits in strArray[4].strip('"').split(",") if splits != ""]
            iNumMultialleles = len(strVars)

            varArray = strArray[9:] # placeholder array of mutations that remains unchanged through iterations

            for i in range (0, (iNumMultialleles)):
                #set ALT allele
                strArray[4] = strVars[i]
                #change variant status
                variant_to_keep = str(i+1)
                # change all number to 0 EXCEPT variant_to_keep
                varArray1 = [variant_to_keep if number == variant_to_keep else '0' for number in varArray.copy()]
                #convert variant_to_keep to 1
                varArray1 = ["1" if number == variant_to_keep else '0' for number in varArray1]
                
                strArray[9:] = varArray1
                fileOutput.write("\t".join(strArray) + "\n")
        else:
            fileOutput.write("\t".join(strArray) + "\n")
print("Done.")

#close the files
fileInput.close()
fileOutput.close()

sys.exit()