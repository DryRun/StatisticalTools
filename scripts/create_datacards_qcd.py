import os
import sys
sys.path.append("/uscms/home/dryu/Dijets/CMSSW_5_3_32_patch3/python/CMSDIJET/QCDAnalysis")
import analysis_configuration_8TeV as analysis_config

datacard_directory = "/uscms/home/dryu/Dijets/data/EightTeeEeVeeBee/Fits/Datacards/condor"
os.chdir(datacard_directory)
#analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM"]
useMCTrigger = False
if useMCTrigger:
	analyses = ["trigbbl_CSVTM", "trigbbh_CSVTM"]
	masses = {"trigbbl_CSVTM":range(350, 850, 50), "trigbbh_CSVTM":range(600, 1250, 50)}
	mjj_min = {"trigbbl_CSVTM":296, "trigbbh_CSVTM":526}
	mjj_max = {"trigbbl_CSVTM":1058, "trigbbh_CSVTM":1607}
else:
	analyses = ["NoTrigger_eta1p7_CSVTM", "NoTrigger_eta2p2_CSVTM"]
	masses = {"NoTrigger_eta1p7_CSVTM":range(350, 850, 50), "NoTrigger_eta2p2_CSVTM":range(600, 1250, 50)}
	mjj_min = {"NoTrigger_eta1p7_CSVTM":296, "NoTrigger_eta2p2_CSVTM":526}
	mjj_max = {"NoTrigger_eta1p7_CSVTM":1058, "NoTrigger_eta2p2_CSVTM":1607}

data_sample = "QCD_TuneZ2star_8TeV_pythia6"

for model in ["Hbb", "RSG"]:
	for analysis in analyses:
		for mass in masses[analysis]:
			data_file_path = analysis_config.get_b_histogram_filename(analysis, data_sample)

			if useMCTrigger:
				signal_pdf_file = analysis_config.get_signal_fit_file(analysis, model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
			else:
				if analysis == "trigbbl_CSVTM" or analysis == "NoTrigger_eta1p7_CSVTM":
					notrig_analysis = "trigbbl_notrig_CSVTM"
				elif analysis == "trigbbh_CSVTM" or analysis == "NoTrigger_eta2p2_CSVTM":
					notrig_analysis = "trigbbh_notrig_CSVTM"
				else:
					print "ERROR : I don't know a no-trigger variant of analysis {}. Please make one, or specify useMCTrigger.".format(analysis) 
					sys.exit(1)
				signal_pdf_file = analysis_config.get_signal_fit_file(notrig_analysis, model, mass, "bukin", interpolated=(not mass in analysis_config.simulation.simulated_masses))
			input_files = [data_file_path, signal_pdf_file]
			command = "python $CMSSW_BASE/src/CMSDIJET/StatisticalTools/scripts/create_datacards_parallel.py {} {}".format(analysis, model)
			command += " --massMin {} --massMax {} --mass {}".format(mjj_min[analysis], mjj_max[analysis], mass)
			command += " --qcd --runFit --condor"
			run_script_path = datacard_directory + "/run_dc_{}_{}_{}.sh".format(model, analysis, mass)
			run_script = open(run_script_path, "w")
			run_script.write("#!/bin/bash\n")
			run_script.write(command + "\n")
			run_script.close()
			condor_command = "csub {} --cmssw --no_retar -F {}".format(run_script_path, ",".join(input_files))
			os.system(condor_command)

			if model == "Hbb" and mass == 750:
				# Run another job with fitBonly, for plotting
				command = "python $CMSSW_BASE/src/CMSDIJET/StatisticalTools/scripts/create_datacards_parallel.py {} {}".format(analysis, model)
				command += " --massMin {} --massMax {} --mass {}".format(mjj_min[analysis], mjj_max[analysis], mass)
				command += " --qcd --runFit --condor --fitBonly"
				run_script_path = datacard_directory + "/run_dc_{}_{}_{}_fitBonly.sh".format(model, analysis, mass)
				run_script = open(run_script_path, "w")
				run_script.write("#!/bin/bash\n")
				run_script.write(command + "\n")
				run_script.close()
				condor_command = "csub {} --cmssw --no_retar -F {}".format(run_script_path, ",".join(input_files))
				os.system(condor_command)



#python create_datacards_parallel.py trigbbh_CSVTM Hbb --massMin 526 --massMax 1607 --massrange 600 1200 50 --fitTrigger --runFit >& log_trigbbh_CSVTM_Hbb_fitTrigger.txt &
#python create_datacards_parallel.py trigbbl_CSVTM Hbb --massMin 296 --massMax 1246 --massrange 350 700 50 --fitTrigger --runFit >& log_trigbbl_CSVTM_Hbb_fitTrigger.txt &
#python create_datacards_parallel.py trigbbh_CSVTM RSG --massMin 526 --massMax 1607 --massrange 600 1200 50 --fitTrigger --runFit >& log_trigbbh_CSVTM_RSG_fitTrigger.txt &
#python create_datacards_parallel.py trigbbl_CSVTM RSG --massMin 296 --massMax 1246 --massrange 350 700 50 --fitTrigger --runFit >& log_trigbbl_CSVTM_RSG_fitTrigger.txt &
