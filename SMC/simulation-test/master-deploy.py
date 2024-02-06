import os
import numpy
import shutil

nrows = 2

for i in range(nrows):
	shutil.copyfile("deploy_template.py", "deploy.py")
	ratio = numpy.random.uniform(low=0.0, high=0.1)
	theta = numpy.random.uniform(low=0.0, high=0.5)
	lamda = theta / ratio

	search_text = "theta          = X"
	replace_text = "theta          = %12.5f" % theta

	search_text_two = "lamda          = X"
	replace_text_two = "lamda          = %12.5f" % lamda

	with open("deploy.py", mode="r") as file:
		data = file.read()
		data = data.replace(search_text, replace_text)

	with open("deploy.py", mode="w") as file:
		file.write(data)

	with open("deploy.py", mode="r") as file:
		datatwo = file.read()
		datatwo = datatwo.replace(search_text_two, replace_text_two)

	with open ("deploy.py", mode="w") as file:
		file.write(datatwo)

	# rename directory
	dirname = "g" + str(i)

	search_text_three = "dirname        = 'g'"
	replace_text_three = "dirname        = '%s'" % dirname

	with open ("deploy.py", mode="r") as file:
		datathree = file.read()
		datathree = datathree.replace(search_text_three, replace_text_three)

	with open ("deploy.py", mode="w") as file:
		file.write(datathree)

	exec(open("./deploy.py").read())
	#os.rename("g", "g"+str(i))

# create master slurm file

f = open("submit-all.sh", "x")
f.write("#!/bin/bash\n")

for i in range(nrows):
	directory = "g"+str(i)
	f.write("cd " + directory + "\n")
	f.write ("sbatch smcslurm.sh\n")
	f.write ("sbatch beastslurm.sh\n")
	f.write ("cd .. \n")

# create master simulate file
f2 = open ("simulate-all.sh", "x")
f2.write("#!/bin/bash\n")

for i in range(nrows):
	directory = "g" + str(i)
	f2.write("cd " + directory + "\n")
	f2.write(". simulate.sh\n")
	f2.write("cd .. \n")

# create master copy data file
f2 = open ("copy-all.sh", "x")
f2.write("#!/bin/bash\n")

for i in range(nrows):
       	directory = "g" + str(i)
       	f2.write("cd " + directory + "\n")
       	f2.write("python3 copydata.py\n")
       	f2.write("cd .. \n")
