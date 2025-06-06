#-------------------------------------------------------------------------------
# Clean the vcf:: drop extra lines and SNPs with incomplete data
# drop variants if MQ is too small or minor allele too rare
# report locations of indels
# Require call from all IM pools
#-------------------------------------------------------------------------------

import sys

vcfile = sys.argv[1] # varscan
samples= 72
MinMQ = 20
MinQual = 20

out2 =open("cleaned."+vcfile, "w") 
out1 =open("indels."+vcfile, "w") 

passfail=[0,0]
src  =open(vcfile, "r") # 
for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t') 
# Chr_04	9888469	.	A	T	999	.	DP=10191;VDB=1.25836e-18;SGB=4574.01;RPB=0.875292;MQB=7.45499e-05;MQSB=0;BQB=0.0262269;MQ0F=0.000294377;ICB=0.00123077;HOB=0.00059453;AC=114;AN=116;DP4=44,16,5007,3197;MQ=57	GT:PL:AD	1/1:255,255,0:11,262	1/1:255,255,0:0,266	1/1:255,255,0:0,238	1/1:255,255,0:1,264	1/1:255,255,0:0,293	1/1:255,255,0:1,284	1/1:255,255,0:0,309	1/1:255,255,0:0,270	1/1:255,255,0:0,289	1/1:255,255,0:0,329	1/1:255,255,0:3,454	1/1:255,255,0:0,498	1/1:255,255,0:1,338	1/1:255,255,0:0,215	1/1:255,255,0:0,392	1/1:255,255,0:2,356	1/1:255,255,0:1,450	1/1:255,255,0:1,538	1/1:255,255,0:2,264	1/1:255,255,0:0,277	1/1:255,255,0:0,611	1/1:255,255,0:2,821	0/1:255,0,198:23,62	0/1:255,0,126:12,311/1:255,45,0:0,15	1/1:46,6,0:0,2	1/1:87,9,0:0,3	./.:0,0,0:0,0	...

	if len(cols)==(samples+9):

		if cols[0]=="#CHROM":
			pass
		elif len(cols[3])==1 and len(cols[4])==1:
			xyt = cols[7].split(";")
			mq = xyt[len(xyt)-1].split("=")
			if mq[0]=="MQ":
				MapQ = int(mq[1])
			else:
				print "fail x",xyt

			if MapQ>=MinMQ and float(cols[5])>=MinQual:

				# require data for all IM samples
				OK=0
				for j in range(9,9+22):
					vv=cols[j].split(":")[0]
					if vv == "./.":
						OK=1
				if OK==0:
					out2.write(line)
		else:
			out1.write(cols[0]+'\t'+cols[1]+'\t'+cols[3]+'\t'+cols[4]+'\n') # Chr_04	9888469	.	A	T


