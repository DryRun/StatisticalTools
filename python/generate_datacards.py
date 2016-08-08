# Take the set of datacards that Tyler produces, and change to different fit functions

models = ["Hbb", "RSG"]
analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM"]
masses = {"trigbbl_CSVTM":range(400, 950, 50), "trigbbh_CSVTM":range(550, 1250, 50)}
datacard_folders = {
	"Hbb":{
		"trigbbl_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidlGcards",
		"trigbbh_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidhGcards",
	},
	"RSG":{
		"trigbbl_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidlRcards",
		"trigbbh_CSVTM":"/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/BiasStudies/davidhRcards",
	},
}
functions = ["f" + str(x) for x in xrange(1, 6)]
parameters = {
	"f1":["p1", "p2", "p3"],
	"f2":["p4", "p5"],
	"f3":["p6", "p7"],
	"f4":["p8", "p9", "p10"],
	"f5":["p11", "p12"]
}
for model in models:
	for analysis in analyses:
		for mass in masses[analysis]:
			for function in functions:
				input_datacard = open(datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + ".txt", 'r')
				output_datacard = open(datacard_folders[model][analysis] + "/datacard_qq_m" + str(mass) + "_" + function + ".txt", "w")
				for line in input_datacard:
					print "[debug] Read in " + line
					processed_line = line
					if "p1" in line or "p2" in line or "p3" in line:
						print "Continuing"
						continue
					if "process" in line and "signal" in line and "background" in line:
						if function != "f1":
							processed_line = processed_line.replace("background", "background" + function.replace("f", ""))
					print "Writing line " + processed_line
					output_datacard.write(processed_line)
				for parameter in parameters[function]:
					output_datacard.write(parameter + "  flatParam\n")
				output_datacard.close()
				input_datacard.close()