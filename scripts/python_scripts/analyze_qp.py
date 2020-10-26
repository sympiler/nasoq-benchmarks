import yaml
import os
import sys
import csv
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as la

header = ["Problem name", "Problem category", "Problem source", "Application", "Number of variables", \
		"QP dimension", "QP norm (order = infty)", "QP condition (2-norm)", "QP min", "QP max", "QP nonzeros", \
			"Number of equality constraints", "Number of inequality constraints", "equality norm (order = infty)", "inequality norm (order = infty)", \
				"equality condition number (2-norm)", "inequality condition number (2-nrom)", "equality nonzeros", "inequality nonzeros", \
					"equality max", "inequality max", "equality min", "inequality min"]

class properties:
	"""
	properties of each Quadratic Problem
	"""
	def __init__(self):
		# the description of problem: name, categorty, application area and source
		self.prob_name = self.prob_type = self.prob_app = self.prob_src = "\\"
		# all norms are evaluated at order infinity
		# and all condition numbers are evaluated at 2-norm

		# number of variables for the Quadratic problem
		self.vars_num = 0
		# the dimension of this Quadratic problem
		self.qp_dim = 0
		# the norm of Quadratic problem matrix
		self.qp_norm = 0
		# the cond number of Quadratic problem matrix
		self.qp_cond = "\\"
		# the number of non-zero entries in Quadratic problem matrix
		self.qp_nnz = 0
		# the max entry of Quadratic problem matrix
		self.qp_max = -np.inf
		# the min entry of Quadratic problem matrix (non-zero)
		self.qp_min = np.inf

		# all data about equality and ineqaulity constraints have default "\\"

		# the number of constraints in equality constraint
		self.equal_constr_num = "\\"
		# the number of constraints in inequality constraint
		self.inequal_constr_num = "\\"
		# the norm of equality constraint matrix
		self.equal_constr_norm = "\\"
		# the norm of inequality constraint matrix
		self.inequal_constr_norm = "\\"
		# the condition number of equality constraint matrix
		self.equal_constr_cond = "\\"
		# the condition number of inequality constraint matrix
		self.inequal_constr_cond = "\\"
		# the number of non-zero entries of equality constraint matrix
		self.equal_constr_nnz = "\\"
		# the number of non-zero entries of inequality constraint matrix
		self.inequal_constr_nnz = "\\"
		# the maximum entry of equality constraint matrix
		self.equal_constr_max = "\\"
		# the maximum entry of inequality constraint matrix
		self.inequal_constr_max = "\\"
		# the minimum entry of equality constraint matrix
		self.equal_constr_min = "\\"
		# the minimum entry of inequality constraint matrix
		self.inequal_constr_min = "\\"
	
	def __str__(self):
		"""
		an analysis report about this problem by the script
		"""
		result = "Problem analyze list:\n"
		result += "	Problem name: {}\n".format(self.prob_name)
		result += "	Problem category: {}\n".format(self.prob_type)
		result += "	Problem source: {}\n".format(self.prob_src)
		result += "	Application: {}\n".format(self.prob_app)
		result += "	Number of variables: " + str(self.vars_num) + "\n"
		result += "	QP dimension: " + str(self.qp_dim) + "\n"
		result += "	QP norm (order = infty): " + str(self.qp_norm) + "\n"
		# result += "	QP condition (2-norm): " + str(self.qp_cond) + "\n"
		result += "	QP min:"  + str(self.qp_min) + "\n"
		result += "	QP max: " + str(self.qp_max) + "\n"
		result += "	QP nonzeros: " + str(self.qp_nnz) + "\n"
		result += "	Number of equality constraints: " + str(self.equal_constr_num) + "\n"
		result += "	Number of inequality constraints: " + str(self.inequal_constr_num) + "\n"
		result += "	equality norm (order = infty): " + str(self.equal_constr_norm) + "\n"
		result += "	inequality norm (order = infty): " + str(self.inequal_constr_norm) + "\n"
		# result += "	equality condition number (2-norm): " + str(self.equal_constr_cond) + "\n"
		# result += "	inequality condition number (2-nrom): " + str(self.inequal_constr_cond) + "\n"
		result += "	equality nonzeros: " + str(self.equal_constr_nnz) + "\n"
		result += "	inequality nonzeros: " + str(self.inequal_constr_nnz) + "\n"
		result += "	equality max: " + str(self.equal_constr_max) + "\n"
		result += "	inequality max: " + str(self.inequal_constr_max) + "\n"
		result += "	equality min: " + str(self.equal_constr_min) + "\n"
		result += "	inequality min: " + str(self.inequal_constr_min) + "\n"

		return result
	
	def map_props(self):
		"""
		make a dictionary to map each property to its label
		"""
		result = {}

		result[header[0]] = self.prob_name
		result[header[1]] = self.prob_type
		result[header[2]] = self.prob_src
		result[header[3]] = self.prob_app
		result[header[4]] = self.vars_num
		result[header[5]] = self.qp_dim
		result[header[6]] = self.qp_norm
		result[header[7]] = self.qp_cond
		result[header[8]] = self.qp_min
		result[header[9]] = self.qp_max
		result[header[10]] = self.qp_nnz
		result[header[11]] = self.equal_constr_num
		result[header[12]] = self.inequal_constr_num
		result[header[13]] = self.equal_constr_norm
		result[header[14]] = self.inequal_constr_norm
		result[header[15]] = self.equal_constr_cond
		result[header[16]] = self.inequal_constr_cond
		result[header[17]] = self.equal_constr_nnz
		result[header[18]] = self.inequal_constr_nnz
		result[header[19]] = self.equal_constr_max
		result[header[20]] = self.inequal_constr_max
		result[header[21]] = self.equal_constr_min
		result[header[22]] = self.inequal_constr_min

		return result

#TODO: I'd suggest rename the following to read_SMP. We can add this to the SMP repo.
def read_SMP(path="test05_0.yml"):
	props = properties()
	with open(path, 'r') as stream:
		data_loaded = yaml.safe_load(stream)
	keys = list(data_loaded.keys())

	# descriptions about the QP
	description = data_loaded[keys[0]]
	description_lst = description.splitlines()
	for feature in description_lst:
		if feature.startswith('category'):
			categorty = feature[feature.find("=") + 2:]
			if categorty:
				props.prob_type = categorty
		elif feature.startswith('application'):
			app = feature[feature.find("=") + 2:]
			if app:
				props.prob_app = app
		elif feature.startswith('name'):
			name = feature[feature.find("=") + 2:]
			if name:
				props.prob_name = name
		elif feature.startswith('source'):
			src = feature[feature.find("=") + 2:]
			if src:
				props.prob_src = src
	# properties of the QP matrix
	qp_str = data_loaded["Quadratic"].splitlines()
	qp_des = qp_str[1]
	qp_des = np.int32(qp_des.split(" "))
	del qp_str[0]; del qp_str[0]

	props.vars_num = qp_des[1]
	props.qp_dim = qp_des[0]
	props.qp_nnz = qp_des[2]

	qp_matrix = sp.lil_matrix((qp_des[0], qp_des[1]))
	# form the QP matrix
	for i in range(qp_des[2]):
		entry = qp_str[i].split(" ")
		qp_matrix[np.int32(entry[0]) - 1, np.int32(entry[1]) - 1] = np.float64(entry[2])

	qp_coo = qp_matrix.tocoo()
	props.qp_norm = la.norm(qp_matrix, ord=np.inf)
	# props.qp_cond = np.linalg.cond(qp_matrix)
	# props.qp_max = np.max(qp_matrix)
	# props.qp_min = np.min(qp_matrix[np.nonzero(qp_matrix)])
	props.qp_max = qp_coo.data[qp_coo.data.argmax()]
	props.qp_min = qp_coo.data[qp_coo.data.argmin()]

	# properties of equality matrix
	if keys.count("Equality"):
		equal_str = data_loaded["Equality"].splitlines()
		equal_des = equal_str[1]
		equal_des = np.int32(equal_des.split(" "))
		del equal_str[0]; del equal_str[0]

		props.equal_constr_num = equal_des[0]
		props.equal_constr_nnz = equal_des[2]

		# form the equality constraint matrix
		equal_matrix = sp.lil_matrix((equal_des[0], equal_des[1]))
		for i in range(equal_des[2]):
			entry = equal_str[i].split(" ")
			equal_matrix[np.int32(entry[0]) - 1, np.int32(entry[1]) - 1] = np.float64(entry[2])
		
		eq_coo = equal_matrix.tocoo()
		props.equal_constr_norm = la.norm(equal_matrix, ord=np.inf)
		# props.equal_constr_cond = np.norm.cond(equal_matrix)
		# props.equal_constr_max = np.max(equal_matrix)
		# props.equal_constr_min = np.min(equal_matrix[np.nonzero(equal_matrix)])
		props.equal_constr_max = eq_coo.data[eq_coo.data.argmax()]
		props.equal_constr_min = eq_coo.data[eq_coo.data.argmin()]
	
	# properties of inequality matrix
	if keys.count("Inequality"):
		inequal_str = data_loaded["Inequality"].splitlines()
		inequal_des = inequal_str[1]
		inequal_des = np.int32(inequal_des.split(" "))
		del inequal_str[0]; del inequal_str[0]

		props.inequal_constr_num = inequal_des[0]
		props.inequal_constr_nnz = inequal_des[2]

		# form the inequality constraint matrix
		inequal_matrix = sp.lil_matrix((inequal_des[0], inequal_des[1]))
		for i in range(inequal_des[2]):
			entry = inequal_str[i].split(" ")
			inequal_matrix[np.int32(entry[0]) - 1, np.int32(entry[1]) - 1] = np.float64(entry[2])
		
		ineq_coo = inequal_matrix.tocoo()
		props.inequal_constr_norm = la.norm(inequal_matrix, ord=np.inf)
		# props.inequal_constr_cond = np.linalg.cond(inequal_matrix)
		# props.inequal_constr_max = np.max(inequal_matrix)
		# props.inequal_constr_min = np.min(inequal_matrix[np.nonzero(inequal_matrix)])
		props.inequal_constr_max = ineq_coo.data[ineq_coo.data.argmax()]
		props.inequal_constr_min = ineq_coo.data[ineq_coo.data.argmin()]

	return props

def main():
	"""
	read all quadratic problems in the repo and
	analyze each of them
	"""
	# parse the variables
	if len(sys.argv) != 3:
		print("analyze_qp.py <path to QPs> <ifMerge>")
		sys.exit(1)
	
	# make the directory for analysis csvs
	csv_folder = sys.argv[1] + "_analysis_tables"
	if not os.path.exists(csv_folder):
		os.makedirs(csv_folder)
	
	# analyze each quadratic problems
	for root, dirs, files in os.walk(sys.argv[1]):

		# setup for reading each file in current directory
		QPs = []
		root_lst = root.split("/")
		if root == ".":
			csv_name = csv_folder + "/smp_repo_QPs_analysis.csv"
		else:
			csv_name = csv_folder + "/{}_QPs_analysis.csv".format(root_lst[-1])

		num_dirs = len(dirs)
		if num_dirs > 0:
			print(str(num_dirs) + " dirs left...")

		# make analysis if the folder is not analyzed before
		if not os.path.exists(csv_name):
			print("analyzing {} ...".format(root))
			num_files = len([file for file in files if (file.endswith(".yml") \
				or file.endswith(".yaml"))])
			#TODO: why do you need to time? We are not interested in this timing. Remove it or explain it to me in slack
			# why you need it.
			start = time.perf_counter()
			for file in files:
				print("{} files left...".format(num_files))
				if file.endswith(".yml") or file.endswith(".yaml"):
					path = root + "/"  + file
					print("reading {} ...".format(path))
					QPs.append(read_SMP(path))
					print("done!")
					num_files -= 1

			# write out a csv for summary of analysis
			if len(QPs) > 0:
				with open(csv_name, "w") as analysis:
					writer = csv.DictWriter(analysis, fieldnames=header)
					writer.writeheader()

					for qp in QPs:
						prop_dict = qp.map_props()
						writer.writerow(prop_dict)
					analysis.flush()
				
				# time the performance for each folder and print
				# adapted from https://stackoverflow.com/questions/27779677/how-to-format-elapsed-time-from-seconds-to-hours-minutes-seconds-and-milliseco
				end = time.perf_counter() - start
				hours, rem = divmod(end, 3600)
				minutes, seconds = divmod(rem, 60)
				print("time consumed: {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))
				print("analyzing {} done!".format(root))

		num_dirs -= 1
#TODO: that would be nice if you have an option to concatenate all csv files into a big one for the whole repo.
	# merge all created csvs into one
	if bool(sys.argv[2]):
		print("begin to merge...")
		csvs = os.listdir(csv_folder)

		# with open(csv_folder + "/merge.csv", "w") as fout:
		# 	i = 0
		# 	while i < len(header):
		# 		fout.write(header[i])
		# 		if i < len(header) - 1:
		# 			fout.write(",")
		# 		i += 1
		# 	# now the rest:    
		# 	for csv_name in csvs:
		# 		f = open(csv_folder + "/" + csv_name, "r")
		# 		# f.nextline() # skip the header
		# 		for line in f:
		# 			fout.write(line)
		# 		f.close()
		combined_csv = pd.concat([pd.read_csv(csv_folder + "/" + f) for f in csvs if os.stat(csv_folder + "/" + f).st_size and f.endswith(".csv")])
		combined_csv.to_csv(csv_folder + "/merged_csvs.csv", index=True)
		print("all csvs are merged")

	return 0

if __name__ == "__main__":
	# props = analyze()
	# print(props)
	main()

# del data_loaded[keys[0]]
# del keys[0]
# with open("{}.txt".format(keys[1]), 'w') as matrixMarket:
# 	matrixMarket.write( data_loaded[keys[1]])

# with open("{}.txt".format(keys[1]), 'r') as matrixMarket:
# 	r = matrixMarket.readlines()
# 	if r[0].startswith("%%"):
# 		del r[0]
# 	s = r[0].splitlines()
# 	s = s[0].split(" ")
# 	result = []

# 	if len(s) == 1:
# 		result = np.array([np.float64(s[0])])
# 	else:
# 		info = []
# 		for id in s:
# 			info.append(np.int32(id))
# 		result = sp.lil_matrix((info[0], info[1]))

# 		del r[0]

# print(s)
# print(result)
# print(result.shape)
# print(info)

# print(len(r))
