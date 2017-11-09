from array import array
from ROOT import TGraph
systematics = {
	"luminosity":0.022,
	"jes":0.01,
	"jer":0.1,
	"bon":0.1,
	"boff":{
		"trigbbl_CSVTM":{
			"Hbb":TGraph(6, array('d',[325,350,400,500,600,750]), array('d', [0.0496984865662,0.0509296944358,0.0527920135972,0.0598252120352,0.0694976071145,0.0920120829753])),
			"ZPrime":TGraph(6, array('d',[325,350,400,500,600,750]), array('d', [0.0495324454323,0.0505803469292,0.0522080272194,0.0589351319081,0.0679690855821,0.0915474509442])),
			"RSG":TGraph(6, array('d',[325,350,400,500,600,750]), array('d', [0.0502374705067,0.0511159167554,0.0531700666596,0.0604638912677,0.0703024641803,0.0954048082064])),
		},
		"trigbbh_CSVTM":{
			"Hbb":TGraph(4, array('d',[600,750,900,1200]), array('d', [0.0700306040447,0.0924340827971,0.115616130354,0.137068438124])),
			"ZPrime":TGraph(4, array('d',[600,750,900,1200]), array('d', [0.0686061711272,0.0919273386348,0.118085204738,0.151517524565])),
			"RSG":TGraph(4, array('d',[600,750,900,1200]), array('d', [0.0707278852876,0.0954554316355,0.119595358442,0.154952753415])),
		}
	},
	"pdfunc":{
		"trigbbl_CSVTM":{
			"Hbb":TGraph(2, array('d', [325,750]), array('d', [0.05, 0.05])),
			"ZPrime":TGraph(2, array('d', [325,750]), array('d', [0.02, 0.03])),
			"RSG":TGraph(2, array('d', [325,750]), array('d', [0.04, 0.04]))
		}, 
		"trigbbh_CSVTM":{
			"Hbb":TGraph(2, array('d', [600, 1200]), array('d', [0.03, 0.03])),
			"ZPrime":TGraph(2, array('d', [600, 1200]), array('d', [0.02, 0.04])),
			"RSG":TGraph(2, array('d', [600, 1200]), array('d', [0.03, 0.03]))
		}
	}
}
