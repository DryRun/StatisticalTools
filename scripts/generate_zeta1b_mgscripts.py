import os
import sys

masses = [325]
masses.extend(range(350, 1250, 50))
basedir = "/home/dryu/Dijets/data/EightTeeEeVeeBee/ZPrime/MG5_aMC_v2_5_5/Zp1b"
final_states = ["u u~", "d d~", "c c~", "s s~", "b b~", "t t~"]
mt = 172.44
#final_states = ["b b~"]
for mass in masses:
	for name, initial_states in {"u":["u u~"], "d":["d d~"], "pp":["p p"]}.iteritems():
		tag = str(mass) + "_" + name
		mg5_script = open("{}/run_{}.txt".format(basedir, tag), "w")
		mg5_script.write("import model Zp -Zp\n")
		first = True
		for initial_state in initial_states:
			for final_state in final_states:
				if final_state == "t t~" and mass < 2 * mt:
					continue
				if first:
					mg5_script.write("generate {} > zph > {}\n".format(initial_state, final_state))
					first = False 
				else:
					mg5_script.write("add process {} > zph > {}\n".format(initial_state, final_state))
		output_directory = "{}/{}".format(basedir, tag)
		mg5_script.write("output {}\n".format(output_directory))
		mg5_script.write("launch\n")
		mg5_script.write("set MZph {}\n".format(mass))
		mg5_script.write("set gz 0.75\n")
		mg5_script.write("set nevents 2500\n")
		mg5_script.write("set ebeam1 4000\n")
		mg5_script.write("set ebeam2 4000\n")
		mg5_script.write("launch {} -i\n".format(output_directory))
		mg5_script.write("print_results --path={}/cross_section.txt --format=short\n".format(output_directory))
		mg5_script.write("exit\n")


		#mg5_script.write("print_results --path={}/results.txt\n".format(output_directory))
		mg5_script.close()
		print "./bin/mg5_aMC {}/run_{}.txt".format(basedir, tag)

