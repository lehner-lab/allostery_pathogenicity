import sys
file=open(sys.argv[1],"r")

for line in file:
	line=line.rstrip("\n")
	line=line.split("\t")

	pos1=line[1]
	pos2=line[0]

	res1_chain=line[2].split("_")[0].split(":")[0]
	res1_atom=line[2].split("_")[0].split(":")[3]
	res2_chain=line[2].split("_")[1].split(":")[0]                       
	res2_atom=line[2].split("_")[1].split(":")[3]	

	ddG_sum=float(line[22])/10
	structure=line[3]
	swapping=line[23]
	
	if structure=="active":
		color="red"
	elif structure=="inactive":
		color="orange"
	elif structure=="non-specific":
		color="blue"

	if swapping=="yes":
		dashes="4"
		
		print("; halfbond = false")
		print("; color = "+color)
		print("; radius = "+str(ddG_sum))
		print("; dashes = "+dashes)
		print("/"+res1_chain+":"+pos1+"@"+res1_atom+" /"+res2_chain+":"+pos2+"@"+res2_atom)

file.close()
