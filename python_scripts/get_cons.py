import pdb,sys,os,urllib
import urllib
pythonversion=sys.version_info[0]

if pythonversion<3:
	import urllib2
	from urllib2 import urlopen
	from urllib import urlencode
else:
	from urllib.request import urlopen
	from urllib.parse import urlencode


def get_conservation(pos):
	url="https://genome.ucsc.edu/cgi-bin/hgTables?"
	form_data={'hgsid':'666516159_sXwdmHvVwfpAowRRzGB5IcIEdQNY',
				'clade':'mammal',
				'org':'human',
				'db':'hg38',
				'hgta_group':'compGeno',
				'hgta_track':'cons100way',
				'hgta_table':'phyloP100way',
				'hgta_regionType':'range',
				'position':pos,
				'hgta_outputType':'wigData',
				'hgta_outFileName':'1',
				'hgta_doTopSubmit':'get output'}
	params=urlencode(form_data).encode('utf-8')
	response=urlopen(url,params)
	data=response.read()
	data=data.decode('utf-8')
	#pdb.set_trace()
	data=data.split('\n')
	data=data[1:]
	out=[]
	#pdb.set_trace()
	for i in data:
		if (len(i)>0) and (i[0]!="#")and (i[0]!='v'):
			ii=i.split()
			ii=[int(ii[0]),float(ii[1])]
			out.append(ii)
	return out


pos="chrY:10000-15000"  # format of genomic range
s=get_conservation(pos) # output  (position,score)

pdb.set_trace()