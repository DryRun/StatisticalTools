# StatisticalTools

Code to perform fits, bias studies, compute limits and significances

## Table of contents

* [Software setup](#software-setup)
* [Limit and significance calculation](#limit-and-significance-calculation)
   * [Resonance shapes](#resonance-shapes)
   * [Datacards](#datacards)
   * [Limit calculation](#limit-calculation)
   * [Significance calculation](#significance-calculation)
   * [Best-fit signal cross section](#best-fit-signal-cross-section)


## Software setup

First, we need to set up out CMSSW working area and check out and compile the `HiggsAnalysis/CombinedLimit` package:

```
setenv SCRAM_ARCH slc6_amd64_gcc491
cmsrel CMSSW_7_4_7
cd CMSSW_7_4_7/src
cmsenv

git clone -b v6.2.0 --depth 1 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

scram b clean
scram b -j8
```

Next, we will also need some dijet specific code and we will put it in the `test/` subdirectory of our CMSSW working area:

```
cd ../test

git clone git://github.com/CMSDIJET/StatisticalTools.git StatisticalTools
```

All of the above steps need to be done only once.

**UPDATE:** As of July 16, 2015 only a private instance of the `StatisticalTools` repository hosted by CERN GitLab is being updated. Here are the steps you should follow to get the latest updates:

```
cd StatisticalTools
git remote add cms-internal ssh://git@gitlab.cern.ch:7999/CMSDIJET/StatisticalTools.git
git fetch cms-internal
git pull --ff-only cms-internal master
```

Alternatively, you can start by directly cloning the private repository

```
git clone ssh://git@gitlab.cern.ch:7999/CMSDIJET/StatisticalTools.git StatisticalTools
```

## Limit and significance calculation

### Resonance shapes

Before we can compute limits or significances, we need signal resonance shapes. Since we'll be using finely binned resonance shapes required by `combine` and RooFit, and given the number of signal mass points, the ROOT files storing resonance shapes will be several MB in size. So in order not to bloat our repositories with many MBs of binary ROOT files, the resonance shapes will not be stored in the `StatisticalTools` repository but will instead be downloaded from the `DijetShapeInterpolator` repository. Nevertheless, it is hard to completely avoid storing binary ROOT files in the repository so in some cases we will still do it (e.g. data dijet spectrum). Generally, this practice should be limited to ROOT files that are small (not more than a few hundred kB) and/or are not expected to change frequently.

To download the resonance shapes, run the following commands:

```
wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/68e8514f4da8849b99b7dfcf1a7834fa55aeefa6/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -P inputs/

wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/68e8514f4da8849b99b7dfcf1a7834fa55aeefa6/ResonanceShapes_qg_13TeV_Scouting_Spring15.root -P inputs/

wget https://github.com/CMSDIJET/DijetShapeInterpolator/raw/68e8514f4da8849b99b7dfcf1a7834fa55aeefa6/ResonanceShapes_qq_13TeV_Scouting_Spring15.root -P inputs/
```

### Datacards

Another essential ingredient for the statistical analysis are datacards and corresponding RooFit workspaces. Here again, we don't necessarily want to store all of these files in the repository since they can be easily remade using scripts available in the repository. Run the following commands to produce datacards for gg, qg, and qq resonances:

```
./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -f gg -o datacards -l 1866 --lumiUnc 0.027 --jesUnc 0.02 --jerUnc 0.1 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2

./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_qg_13TeV_Scouting_Spring15.root -f qg -o datacards -l 1866 --lumiUnc 0.027 --jesUnc 0.02 --jerUnc 0.1 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2

./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_qq_13TeV_Scouting_Spring15.root -f qq -o datacards -l 1866 --lumiUnc 0.027 --jesUnc 0.02 --jerUnc 0.1 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2
```

For more command-line options, run

```
./scripts/createDatacards.py -h
```

For Bayesian limits a separate set of datacards and RooFit workspaces is needed. This is because the background function parameters either need to be fixed (the `--fixBkg` option), if ignoring the background systematics, or need to be decorrelated (the `--decoBkg` option), if taking the background systematics into account:

```
./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_gg_13TeV_Scouting_Spring15.root -f gg -o datacards_Bayesian -l 1866 --lumiUnc 0.027 --jesUnc 0.02 --jerUnc 0.1 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2 --decoBkg

./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_qg_13TeV_Scouting_Spring15.root -f qg -o datacards_Bayesian -l 1866 --lumiUnc 0.027 --jesUnc 0.02 --jerUnc 0.1 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2 --decoBkg

./scripts/createDatacards.py --inputData inputs/rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root --dataHistname mjj_mjjcor_gev --inputSig inputs/ResonanceShapes_qq_13TeV_Scouting_Spring15.root -f qq -o datacards_Bayesian -l 1866 --lumiUnc 0.027 --jesUnc 0.02 --jerUnc 0.1 --massrange 1000 1500 50 --runFit --p1 5 --p2 7 --p3 0.4 --massMin 838 --massMax 2037 --fitStrategy 2 --decoBkg
```

### Limit calculation

We can now run `combine` to compute the limits. In the following examples the Asymptotic CL<sub>S</sub> method is used:

```
./scripts/runCombine.py -M Asymptotic -d datacards -f gg --massrange 1000 1500 50

./scripts/runCombine.py -M Asymptotic -d datacards -f qg --massrange 1000 1500 50

./scripts/runCombine.py -M Asymptotic -d datacards -f qq --massrange 1000 1500 50
```

For more command-line options, run

```
./scripts/runCombine.py -h
```

In particular, note the `--condor` option which allows you to process all mass points in parallel using the Condor batch system (currently only Condor at FNAL is supported).

To produce the final limit plots, run

```
./scripts/plotLimits.py -M Asymptotic -l logs -f gg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotLimits.py -M Asymptotic -l logs -f qg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotLimits.py -M Asymptotic -l logs -f qq --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults
```

For more command-line options, run

```
./scripts/plotLimits.py -h
```

If you are only interested in producing or modifying plots using already computed results, you can use the `-r` (`--results_file`) instead of the `-l` (`--logs_path`) option as in the following example:

```
./scripts/plotLimits.py -M Asymptotic -r results/limits_Asymptotic_gg_Run2_Scouting_13TeV_DATA_1p9_invfb.py -f gg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps
```

To repeat the same procedure for Bayesian limits (only observed limits for now), run the following set of commands

```
./scripts/runCombine.py -M MarkovChainMC -d datacards_Bayesian -f gg --massrange 1000 1500 50

./scripts/runCombine.py -M MarkovChainMC -d datacards_Bayesian -f qg --massrange 1000 1500 50

./scripts/runCombine.py -M MarkovChainMC -d datacards_Bayesian -f qq --massrange 1000 1500 50


./scripts/plotLimits.py -M MarkovChainMC -l logs -f gg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotLimits.py -M MarkovChainMC -l logs -f qg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotLimits.py -M MarkovChainMC -l logs -f qq --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults
```

### Significance calculation

For the significance calculation we again use the `runCombine.py` script as in the following examples (note the `--signif` option):

```
./scripts/runCombine.py -M ProfileLikelihood --signif -d datacards -f gg --massrange 1000 1500 50

./scripts/runCombine.py -M ProfileLikelihood --signif -d datacards -f qg --massrange 1000 1500 50

./scripts/runCombine.py -M ProfileLikelihood --signif -d datacards -f qq --massrange 1000 1500 50
```

To produce the final significance plots, run

```
./scripts/plotSignificance.py -M ProfileLikelihood -l logs -f gg --massrange 1000 1500 50 --sigRange 2.5 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotSignificance.py -M ProfileLikelihood -l logs -f qg --massrange 1000 1500 50 --sigRange 2.5 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotSignificance.py -M ProfileLikelihood -l logs -f qq --massrange 1000 1500 50 --sigRange 2.5 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults
```

For more command-line options, run

```
./scripts/plotSignificance.py -h
```

If you are only interested in producing or modifying plots using already computed results, you can use the `-r` (`--results_file`) instead of the `-l` (`--logs_path`) option as in the following example:

```
./scripts/plotSignificance.py -M ProfileLikelihood -r results/significance_ProfileLikelihood_gg_Run2_Scouting_13TeV_DATA_1p9_invfb.py -f gg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps
```

### Best-fit signal cross section

The `runCombine.py` script can also be used to extract the best-fit signal cross sections as in the following examples:

```
./scripts/runCombine.py -M MaxLikelihoodFit -d datacards -f gg --rMin -50 --rMax 50 --massrange 1000 1500 50

./scripts/runCombine.py -M MaxLikelihoodFit -d datacards -f qg --rMin -50 --rMax 50 --massrange 1000 1500 50

./scripts/runCombine.py -M MaxLikelihoodFit -d datacards -f qq --rMin -50 --rMax 50 --massrange 1000 1500 50
```

To produce the final cross section plots, run

```
./scripts/plotSignalXSec.py -l logs -f gg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotSignalXSec.py -l logs -f qg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults

./scripts/plotSignalXSec.py -l logs -f qq --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps --printResults
```

For more command-line options, run

```
./scripts/plotSignalXSec.py -h
```

If you are only interested in producing or modifying plots using already computed results, you can use the `-r` (`--results_file`) instead of the `-l` (`--logs_path`) option as in the following example:

```
./scripts/plotSignalXSec.py -r results/signal_xs_MaxLikelihoodFit_gg_Run2_Scouting_13TeV_DATA_1p9_invfb.py -f gg --massrange 1000 1500 50 --extraText Preliminary --lumi_sqrtS="1.9 fb^{-1} (13 TeV)" --postfix Run2_Scouting_13TeV_DATA_1p9_invfb --fileFormat eps
```
