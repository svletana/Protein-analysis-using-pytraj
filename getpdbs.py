import pypdb
import os

pdbs = ['1ubq']
for p in pdbs:
	try:
		pdb_file = pypdb.get_pdb_file(p)
		f = open('{}.pdb'.format(p), 'w')
		f.write(pdb_file)
		f.close()
		os.system("mkdir {}".format(p))
		os.system("mkdir {}/input".format(p))
		os.system("mv {}.pdb {}/input".format(p, p))
		print("{} DONE".format(p))
	except:
		print('Failed to get pdb file')
