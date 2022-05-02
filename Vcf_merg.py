# Pieriandx Assignment I
# VCF files merging script using Python 2.7
# Author: Shyam Ingle
# Date: 09/01/2022

###################################################################

#Finding and merging common calls and unique calls form both files

def split_vcf_file(vcfdata):
	#Get VCF call by spliting vcf file 
	split_data = vcfdata.split("#CHROM")
	#Return comments and calls array
	return split_data
	

def get_calls(call_data):
	#Get calls array by split newline
	calls = call_data.split("\n")
	return calls[1:len(calls)-1]


def merg_and_get_uniq_calls(call_lines):
	#Merg chr with position and calls (Merg CHROM+ POS + REF + ALT)
	uniq_calls = []
	for line in call_lines: 
		line_array = line.split("\t")
		uniq_calls.append(line_array[0]+"_"+line_array[1]+"_"+line_array[3]+"_"+line_array[4])
	return uniq_calls


def get_matched_line(calls,lines_array): 
    # Find common call lines from both VCF files
    ca = calls.split("_")
    for line in lines_array:
        la = line.split("\t")
        if ca[0] == la[0]:
            if ca[1] == la[1]:
                if ca[2] == la[3]:
                    if ca[3] == la[4]:
			            #log.write("Common calls: "+ str(line.split("\t")[:6])+"\n") # print common calls in log file
                        return line
                    else:
                        return "NULL"


###################################################################

import sys
import argparse


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description='VCF files Merging Programm')
	parser.add_argument("-v1", "--VCF1", help="First VCF file name is required", required=True, action="store")
	parser.add_argument("-v2", "--VCF2", help="Second VCF file name is required", required=True, action="store")
	parser.add_argument("-l", "--log", help="Log file name need to mention", required=True, action="store")
	args = parser.parse_args()



log = open(sys.argv[6],'w')# log file 

if sys.argv[2].endswith('.vcf'):
	file1 = open(sys.argv[2],'r') # First vcf file argument
	vcf1 = file1.read() # Reading VCF first file
	log.write(sys.argv[2]+" is opened sucessfully!\n")
	file1.close()
else:
	print "-v1 File is not VCF type!"
	log.write("-v1 File is not VCF type!\n")
	exit()
	
if sys.argv[4].endswith('.vcf'):
	file2 = open(sys.argv[4],'r')# Second cvf file argument
	vcf2 = file2.read() # Reading VCF second file
	log.write(sys.argv[4]+" is opened sucessfully!\n")
	file2.close()
else:
	print "-v2 File is not VCF type!"
	log.write("-v2 File is not VCF type!\n")
	exit()


#Getting uniq calls from vcf1 
data1 = split_vcf_file(vcf1)
call_lines1 = get_calls(data1[1])
uniq_calls1 = merg_and_get_uniq_calls(call_lines1)


#Getting uniq calls from vcf2 
data2 = split_vcf_file(vcf2)
call_lines2 = get_calls(data2[1])
uniq_calls2 = merg_and_get_uniq_calls(call_lines2)




###############################################################

#Finding common tags and renamed with regarding software name

#Getting headers from vcf1 and vcf2

headers1 = data1[0].split("\n")
headers2 = data2[0].split("\n")


#Getting INFO and FORMAT tags from vcf1
INFO_tags1 = []; FORMAT_tags1 = []
for hline in headers1:
	if '##INFO=<ID=' in hline:
		hline_array = hline.split(',')
		INFO_tags1.append(hline_array[0])

	if '##FORMAT=<ID=' in hline:
		hline_array = hline.split(',')
		FORMAT_tags1.append(hline_array[0])


#Getting INFO and FORMAT tags from vcf2
INFO_tags2 = []; FORMAT_tags2 = []
for hline in headers2:
	if '##INFO=<ID=' in hline:
		hline_array = hline.split(',')
		INFO_tags2.append(hline_array[0])
		#print hline_array[0]

	if '##FORMAT=<ID=' in hline:
		hline_array = hline.split(',')
		FORMAT_tags2.append(hline_array[0])
		#print hline_array[0]

		
#Finding common tags from INFO and FORMAT

common_INFO =  list(set(INFO_tags1).intersection(INFO_tags2))

if len(common_INFO) != 0:
	for info_tag in common_INFO:
		head1 = data1[0].replace(info_tag,info_tag[:11]+'Freebayes_'+info_tag[11:])
		head2 = data2[0].replace(info_tag,info_tag[:11]+'VarScan_'+info_tag[11:])
	log.write("Common headers INFO tags found in both VCF files.\n")

else:
	head1 = data1[0]
	head2 = data2[0]
	log.write("No Common headers INFO tags found in both VCF files.\n")



common_FORMAT =  list(set(FORMAT_tags1).intersection(FORMAT_tags2))

if len(common_FORMAT) != 0:
	for format_tag in common_FORMAT:
		head1 = head1.replace(format_tag,format_tag[:13]+'Freebayes_'+format_tag[13:])
		head2 = head2.replace(format_tag,format_tag[:13]+'VarScan_'+format_tag[13:])
	log.write("Common headers FORMAT tags found in both VCF files.\n")

else:
	head1 = data1[0]
	head2 = data2[0]
	log.write("No Common headers FORMAT tags found in both VCF files.\n")


#Merging renamed tags Headers from both VCF Files
new_header = head1+head2
log.write("Headers are merged from both VCF files.\n")
#################################################################


#Get renamed INFO and FORMAT tag calls from vcf1  
newvcf1calls = []
for line in call_lines1:
    line_array = line.split("\t")
    uqc1 = (line_array[0]+"_"+line_array[1]+"_"+line_array[3]+"_"+line_array[4])
    if uqc1 in uniq_calls2:
        matched_line =  get_matched_line(uqc1,call_lines2)
        calls2_line_array = matched_line.split("\t")
        line_array[7] = line_array[7]+calls2_line_array[7]
        #print line_array[7] 
        line_array[7] = line_array[7]+';calledBy=Freebayes+VarScan'
        for format_tag in common_FORMAT:		
            line_array[8] = line_array[8].replace(format_tag[13:],'Freebayes_'+format_tag[13:])
            calls2_line_array[8] = calls2_line_array[8].replace(format_tag[13:],'VarScan_'+format_tag[13:])
        line_array[8] = line_array[8]+calls2_line_array[8]
        newvcf1calls.append(line_array)
	#log.write("Common calls Identifed and tags are renemed.\n")
    else:
		line_array[7] = line_array[7]+';calledBy=Freebayes'
		newvcf1calls.append(line_array)


#Get renamed INFO & FORMAT tags and uniq calls from vcf2  
newvcf2calls = []
for line in call_lines2:
	line_array = line.split("\t")
    if not (line_array[0]+"_"+line_array[1]+"_"+line_array[3]+"_"+line_array[4]) in uniq_calls1:
    line_array[7] = line_array[7]+';calledBy=VarScan'
	newvcf2calls.append(line_array)
log.write("Common calls Identifed and tags are renamed.\n")
		

#Merging both files calls

for uniq_vcf2_calls in  newvcf2calls:
	newvcf1calls.append(uniq_vcf2_calls)
   
###########################################################################   
   
#Writing Creating Merg VCF file

new_calls = ""

for line in newvcf1calls:
    k = ""
    for wc in line: k=k+wc+"\t"
    new_calls = (new_calls+k+"\n")

OUT = open("merged_"+sys.argv[2][:len(sys.argv[2])-4]+"_"+sys.argv[4][:len(sys.argv[4])-4]+".vcf",'w')
OUT.write(new_header+"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tunknown\n"+new_calls)
print "**Merged VCF file Created successfully!**"
log.write("**Merged VCF file Created successfully!**")
log.close()
OUT.close()

###########################################################################
























