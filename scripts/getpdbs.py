import pypdb
import os

pdbs = ['1ubq']
for p in pdbs:
	pdb_file = pypdb.get_pdb_file(p)
	with open('{}.pdb'.format(p), 'w') as f:
		f.write(pdb_file)
	os.system("mkdir {}".format(p))
	os.system("mkdir {}/input".format(p))
	os.system("mv {}.pdb {}/input".format(p, p))
	print("{} done".format(p))
