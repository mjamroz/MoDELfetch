#!/usr/bin/env python
from re import compile
import urllib
# wget "http://mmb.pcb.ub.es/MoDEL/index.jsf?=&idSimulation=189"
# public domain MoDEL parser

p = compile(r'<.*?>')
def remove_html_tags(data):
    global p
    return p.sub('', data.strip())

def fetchData(idx):
	show=False
	showit=0
	f = urllib.urlopen("http://mmb.pcb.ub.es/MoDEL/index.jsf?=&idSimulation="+str(idx))
	s = f.read()
	f.close()

	data = []
	d = s.split("\n")
	jestTraj=False
	for line in d:
		if "type=ca" in line:
			jestTraj=True
			break
	if not jestTraj:
		return []

	for line in d:

		if show and len(line)>2:
			if "</td>" in line or '</label>' in line:
				stripped = remove_html_tags(line)
				if stripped!="":
					showit+=1
					data.append(stripped)
				
		if "<td>PDB code</td>" in line:
			show=True
		if showit==11:
			break
			

	for i in range(len(d)):
		if "<td>Simulation time</td>" in d[i]:
			data.append("time")
			data.append(remove_html_tags(d[i+1]))
			break
	
	data.append("http://mmb.pcb.ub.es/MoDEL/servlet/trajectory?idSimulation="+str(idx)+"&type=ca")
	if(len(data)<3):data=[]
	return data 
#['1AK2', 'Program Version', 'AMBER 8.0', 'Force field', 'parm99 (standard amber forcefield)', 'Salt Concentration', '52.99 mM', 'Total atoms', '37728', 'Total residues', '220', 'Simulation time', '10000 ps', 'http://mmb.pcb.ub.es/MoDEL/servlet/trajectory?idSimulation=230&type=ca']

fw = open("MoDEL.data.txt","w")
fw.write("#%3s %10s %10s %5s %80s # %s\n" %("pdb","software","time","len","trajectory","ff"))

for i in range(80,5000):
	d = fetchData(i)
	if (len(d)>1): 
		#simlen = int(d[12].split()[0])
		
		#if (simlen==10000): # only 10ns
#		print "%4s %10s %10s %5s %80s # %s\n" %(d[0],d[2],d[12],d[10],d[-1],d[4])
		fw.write("%4s %10s %10s %5s %80s # %s\n" %(d[0],d[2],d[12],d[10],d[-1],d[4]))
fw.close()
