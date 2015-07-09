# StatisticalTools

Code to perform fits, bias studies, compute limits and significances

## Table of contents

* [Software setup](#software-setup)
* [Limit and significance calculation](#limit-and-significance-calculation)
   * [Resonance shapes](#resonance-shapes)
   * [Datacards](#datacards)
   * [Limit calculation](#limit-calculation)
   * [Significance calculation](#significance-calculation)


## Software setup

First, we need to set up out CMSSW working area and check out and compile the `HiggsAnalysis/CombinedLimit` package:

```
setenv SCRAM_ARCH slc6_amd64_gcc481
cmsrel CMSSW_7_1_5
cd CMSSW_7_1_5/src
cmsenv

git clone -b v5.0.1 git://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit

scram b
```

Next, we will also need some dijet specific code and we will put it in the `test/` subdirectory of our CMSSW working area:

```
cd ../test

git clone git://github.com/CMSDIJET/DijetShapeInterpolator.git DijetShapeInterpolator
git clone git://github.com/CMSDIJET/StatisticalTools.git StatisticalTools

```

All of the above steps need to be done only once.

## Limit and significance calculation

### Resonance shapes

Before we can compute limits or significances, we need signal resonance shapes. Since we'll be using finely binned resonance shapes required by `combine` and RooFit and given the number of signal mass points, the ROOT files storing resonance shapes will be several MB in size. So in order not to bloat our repositories with many MBs of binary ROOT files, the resonance shapes will not be stored in the `DijetShapeInterpolator` or `StatisticalTools` repositories but will instead be produced using scripts and information stored in the `DijetShapeInterpolator` repository. Nevertheless, it is hard to completely avoid storing binary ROOT files in the repository so in some cases we will still do it (e.g. data dijet spectrum). Generally, this practice should be limited to ROOT file that are small (not more than a few hundred kB) and/or are not expected to change frequently.

To produce the resonance shapes, go to the `DijetShapeInterpolator` package

```
cd DijetShapeInterpolator
```

and run the following commands:

```
./getResonanceShapes.py -i inputs/input_shapes_gg_13TeV_PU30_Spring15.py -f gg --fineBinning --massrange 500 9000 100 -o ResonanceShapes_gg_13TeV_PU30_Spring15.root
./getResonanceShapes.py -i inputs/input_shapes_qg_13TeV_PU30_Spring15.py -f qg --fineBinning --massrange 500 9000 100 -o ResonanceShapes_qg_13TeV_PU30_Spring15.root
./getResonanceShapes.py -i inputs/input_shapes_qq_13TeV_Phys14Spring15Mix.py -f qq --fineBinning --massrange 500 9000 100 -o ResonanceShapes_qq_13TeV_Phys14Spring15Mix.root
```

This will produce gg, qg, and qq resonance shapes. For more command line options, run

```
./getResonanceShapes.py -h
```

Finally, move the produced resonance shapes to the `StatisticalTools` package:

```
mv ResonanceShapes_gg_13TeV_PU30_Spring15.root ResonanceShapes_qg_13TeV_PU30_Spring15.root ResonanceShapes_qq_13TeV_Phys14Spring15Mix.root ../StatisticalTools/inputs/
cd ../StatisticalTools
```

### Datacards

Another essential ingredient for the statistical analysis are datacards and corresponding RooFit workspaces. Here again, we don't necessarily want to store all of these files in the repository since they can be easily remade using scripts available in the repository. Run the following commands to produce datacards for gg, qg, and qq resonances:

```
./scripts/createDatacards.py --inputData inputs/dijetFitResults_FuncType0_nParFit4_MC_1invfb.root --dataHistname hist_mass_1GeV --inputSig inputs/ResonanceShapes_gg_13TeV_PU30_Spring15.root -f gg -o datacards -l 1000 --massrange 1300 7000 100 --runFit
./scripts/createDatacards.py --inputData inputs/dijetFitResults_FuncType0_nParFit4_MC_1invfb.root --dataHistname hist_mass_1GeV --inputSig inputs/ResonanceShapes_qg_13TeV_PU30_Spring15.root -f qg -o datacards -l 1000 --massrange 1300 7000 100 --runFit
./scripts/createDatacards.py --inputData inputs/dijetFitResults_FuncType0_nParFit4_MC_1invfb.root --dataHistname hist_mass_1GeV --inputSig inputs/ResonanceShapes_qq_13TeV_Phys14Spring15Mix.root -f qq -o datacards -l 1000 --massrange 1300 7000 100 --runFit
```

For more command line options, run

```
./scripts/createDatacards.py -h
```

### Limit calculation

We can now run `combine` to compute the limits. In the following examples the Asymptotic CL<sub>S</sub> method is used:

```
./scripts/runCombine.py -M Asymptotic -d datacards -f gg --massrange 1300 7000 100
./scripts/runCombine.py -M Asymptotic -d datacards -f qg --massrange 1300 7000 100
./scripts/runCombine.py -M Asymptotic -d datacards -f qq --massrange 1300 7000 100
```

For more command line options, run

```
./scripts/runCombine.py -h
```

To produce the final limit plots, run:

```
./scripts/plotLimits.py -l logs -f gg --massrange 1300 7000 100
./scripts/plotLimits.py -l logs -f qg --massrange 1300 7000 100
./scripts/plotLimits.py -l logs -f qq --massrange 1300 7000 100
```

For more command line options, run

```
./scripts/plotLimits.py -h
```

If you are only interested in producing or modifying plots using already computed results, you can use the `-r` (`--results_file`) instead of the `-l` (`--logs_path`) option as in the following example:

```
./scripts/plotLimits.py -r results/limits_qq_Run2_13TeV_MC.py -f qq --massrange 1300 7000 100
```

### Significance calculation

For the significance calculation we again use the `runCombine.py` script as in the following examples (note the `--signif` option):

```
./scripts/runCombine.py -M ProfileLikelihood --signif -d datacards -f gg --massrange 1300 6000 100
./scripts/runCombine.py -M ProfileLikelihood --signif -d datacards -f qg --massrange 1300 6000 100
./scripts/runCombine.py -M ProfileLikelihood --signif -d datacards -f qq --massrange 1300 6000 100
```

To produce the final significance plots, run:

```
./scripts/plotSignificance.py -l logs -f gg --massrange 1300 6000 100
./scripts/plotSignificance.py -l logs -f qg --massrange 1300 6000 100
./scripts/plotSignificance.py -l logs -f qq --massrange 1300 6000 100
```

For more command line options, run

```
./scripts/plotSignificance.py -h
```

If you are only interested in producing or modifying plots using already computed results, you can use the `-r` (`--results_file`) instead of the `-l` (`--logs_path`) option as in the following example:

```
./scripts/plotSignificance.py -r results/significance_qq_Run2_13TeV_MC.py -f qq --massrange 1300 6000 100
```
