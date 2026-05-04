import os

indir = "ValidationInstances"
outdir = "Parsed_Large_Instances"
file_list = os.listdir(indir)
filename = file_list[0]

def parse_instance(filename):
	filepath = indir + '/' + filename
	infile = open(filepath, 'r')
	outfile = open(outdir + '/' + filename + "_parsed", 'a')
	tmp = infile.read()
	tmp = tmp.split('\n')
	count = 0
	products = 0
	for l in tmp:
		if 'Periods' in l:
			buff = l.split()
			periode = int(buff[-1][:-1])
			outfile.write(str(periode) + "\n")
			# print(periode)

		if 'Items' in l:
			buff = l.split()
			products = int(buff[-1][:-1])
			outfile.write(str(products) + '\n')
			# print(products)

		if '|' in l and count<products:
			count += 1
			buff = l.split(', ')
			buff[0] = buff[0].split('|')[-1]
			file_buff = ''
			for c in buff:
				file_buff += str(c) + " "
			outfile.write(file_buff + '\n')
			# print(file_buff)

for filename in file_list:
	parse_instance(filename)