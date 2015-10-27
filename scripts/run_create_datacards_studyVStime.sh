outdir_datacard="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_5/src/StatisticalTools/datacards_dataRunD_studyVStime/" 
outdir="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_1_5/src/StatisticalTools/scripts/results_fitSB_RunD_studyVStime/" 

name="DCSonly_sample78_fitTo6500"
#name="golden_fullPeriod"
#name="golden_period1"
#name="golden_period2"
#name="golden_period3"
#name="golden_periodA"
#name="golden_periodB"
#name="golden_periodC"
#name="DCSonly_fullPeriod"
#name="DCSonly_period1"
#name="DCSonly_period2"
#name="DCSonly_period3"

lumi=789
#lumi=547
xsec=0.8333E-01
inputsig="dcap://cmsrm-se01.roma1.infn.it//pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/mc/Spring15_JEC_Summer15_25nsV5_testsignal_20151016_104616/rootfile_QstarToJJ_M_4000_TuneCUETP8M1_13TeV_pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1__MINIAODSIM_JEC_Summer15_25nsV5_20151016_104616_0_reduced_skim.root"
#inputbkg="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#inputbkg="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period1//histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period1//histo_data_mjj_fromTree.root"
#inputbkg="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period2//histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period2//histo_data_mjj_fromTree.root"
#inputbkg="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period3//histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period3//histo_data_mjj_fromTree.root"
#inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF/histo_data_mjj_fromTree.root"
#inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period1/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period1/histo_data_mjj_fromTree.root"
#inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period2/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period2/histo_data_mjj_fromTree.root"
#inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period3/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period3/histo_data_mjj_fromTree.root"
#inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodA/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodA/histo_data_mjj_fromTree.root"
#inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodB/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodB/histo_data_mjj_fromTree.root"
#inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodC/histo_data_mjj_fromTree.root"
#inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodC/histo_data_mjj_fromTree.root"
inputbkg=" /cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_JEC_Summer15_25nsV5_withSF_sample78/histo_data_mjj_fromTree.root"
inputdata="/cmshome/gdimperi/Dijet/CMSDIJETrepo/CMSSW_7_4_3_Dijet/src/CMSDIJET/DijetRootTreeAnalyzer/scripts/plots_data4T_Run2015D_DCSonly_JEC_Summer15_25nsV5_withSF_sample78/histo_data_mjj_fromTree.root"

cmd_datacard="python create_datacard.py --histpdfSig --fitDat -m 4000 --lumi $lumi -n $name --sigEff 1. --sigXS $xsec -o $outdir_datacard --inputSig $inputsig --inputData $inputdata --inputBkg $inputbkg"
cmd_combine="combine -M MaxLikelihoodFit $outdir_datacard/Qstar4000_datacard_$name.txt --rMin -10000 --rMax 100 --saveNormalizations --saveWorkspace --verbose=2"
cmd_mv="mv mlfit.root  $outdir/mlfit_$name.root"
cmd_plot="python plotFit_data.py --inputFileData $inputdata --inputFileRes $outdir/mlfit_$name.root  --inputWS $outdir_datacard/Qstar4000_workspace_$name.root   --lumi $lumi --xsec $lumi -o  $outdir --tag $name"

#echo $cmd_datacard
#echo $cmd_combine
#echo $cmd_mv
#echo $cmd_plot

mkdir -p $outdir_datacard
mkdir -p $outdir

$cmd_datacard
echo ""
echo "****************************** combine output **************************************"
$cmd_combine
echo "**************************** end combine output ************************************"
echo ""
$cmd_mv
$cmd_plot

#
#
#
#
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 12:57 plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period2
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:04 plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period1
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:04 plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_period3
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:06 plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodA
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:09 plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodB
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:09 plots_data4T_Run2015D_golden_547pb-1_JEC_Summer15_25nsV5_withSF_periodC
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:12 plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period1
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:14 plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period2
#drwxr-xr-x 2 gdimperi cms 4.0K Oct 17 13:14 plots_data4T_Run2015D_DCSonly_974pb-1_JEC_Summer15_25nsV5_withSF_period3

