f = open("resultROADEF3.txt")
buff = f.read()
buff = buff.split("\n")
buff = [l.split() for l in buff]

gammaval = [1,21,41,61,81]
tauval = [0, 2, 4, 6, 8, 10]

resK = {}
resS = {}

for l in buff[:-1]:
	if l[0] == 'KC' and int(l[1]) in gammaval and int(l[2]) in tauval:
		a = '{0:.3}'.format(float(l[3]))
		b = '{0:.3}'.format(float(l[4])/float(l[1]))
		print("KC", l[1], l[2], a + " | " + b)
		resK[(l[1], l[2])] = (a,b)

for l in buff[:-1]:
	if l[0] == 'STANDARD' and int(l[1]) in gammaval and int(l[2]) in tauval:
		a = '{0:.3}'.format(float(l[3]))
		b = '{0:.3}'.format(float(l[4])/float(l[1]))
		print("STANDARD",l[1], l[2], a + " | " + b)
		resS[(l[1], l[2])] = (a,b)

for g in gammaval:
	for t in tauval:
		print (str(g),str(t)), resK[(str(g),str(t))], "KC - STANDARD", resS[(str(g),str(t))]