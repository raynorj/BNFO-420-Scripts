import sys

def LoadCSV(filepath, sep=','):
	fh = open(filepath, "r")
	data = fh.readlines()
	fh.close()

	rows = []
	for i in range(0, len(data)):
		rows.append(data[i].replace("\"", "").split(sep))

	for i in range(len(rows)):
		for j in range(len(rows[i])):
			rows[i][j] = rows[i][j].strip()
		
	return rows

m24 = LoadCSV(sys.argv[1], '\t')

#start forming R command
command = "graph = make_graph(~"

#iterate through each gene
for i in range(1, len(m24)):

	#begin argument with current gene
	gene = "\""+m24[i][0]+"\""+"-"
	argument = "\""+m24[i][0]+"\""+"-"
	
	#for each other gene, verify that the p-value is below the significance value 
	#if so, add that gene to the current argument
	for j in range(1, len(m24[i])):
		if (float(m24[i][j]) < 0.01) and (i != j):
			argument += "\""+m24[0][j]+"\"" + ":"

	#check to see if no neighbors exist	
	if (argument == gene):
		continue

	#end argument with comma
	argument = argument[0:-1] + ",\n"

	command += argument
	
#end command with end parenthesis
command = command[0:-2] + ")"

print(command)