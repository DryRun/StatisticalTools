import os
import sys

#masses = [1250., 1500., 1750., 2000., 2500., 3000., 4000.]
masses = [1750.]
basedir = "/home/dryu/Dijets/data/EightTeeEeVeeBee/ZPrime/MG5_aMC_v2_5_5/ZpATLAS"
os.system("mkdir -pv " + basedir)
final_states = ["u u~", "d d~", "c c~", "s s~", "b b~", "t t~"]
mt = 172.44
#final_states = ["b b~"]
for mass in masses:
	for name, initial_states in {"u":["u u~"], "d":["d d~"], "pp":["p p"]}.iteritems():
		tag = str(mass) + "_" + name
		mg5_script = open("{}/run_{}.txt".format(basedir, tag), "w")
		mg5_script.write("import model Zp2 -Zp2\n")
		first = True
		for initial_state in initial_states:
			for final_state in final_states:
				if final_state == "t t~" and mass < 2 * mt:
					continue
				if first:
					mg5_script.write("generate {} > zp > {}\n".format(initial_state, final_state))
					first = False 
				else:
					mg5_script.write("add process {} > zp > {}\n".format(initial_state, final_state))
		output_directory = "{}/{}".format(basedir, tag)
		mg5_script.write("output {}\n".format(output_directory))
		mg5_script.write("launch\n")
		mg5_script.write("set MZP {}\n".format(mass))
		mg5_script.write("set nevents 2500\n")
		mg5_script.write("set ebeam1 6500\n")
		mg5_script.write("set ebeam2 6500\n")
		mg5_script.write("launch {} -i\n".format(output_directory))
		mg5_script.write("print_results --path={}/cross_section.txt --format=short\n".format(output_directory))
		mg5_script.write("exit\n")


		#mg5_script.write("print_results --path={}/results.txt\n".format(output_directory))
		mg5_script.close()
		print "./bin/mg5_aMC {}/run_{}.txt".format(basedir, tag)

