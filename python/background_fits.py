import CMSDIJET.StatisticalTools.trigger_efficiency as trigger_efficiency

def BackgroundFit_f1(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3))))

def BackgroundFit_f2(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2]

def BackgroundFit_f3(x, par):
	return par[0] / (1 + par[1] * (x[0] / 8.e3))**par[2]

def BackgroundFit_f4(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 <= 0:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) < 1.e-15:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3])

def BackgroundFit_f5(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2]

def BackgroundFit_f6(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 + par[3] * (x[0]/8.e3)**3 <= 0:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2  + par[3] * (x[0]/8.e3)**3)**par[4]) < 1.e-15:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 + par[3] * (x[0]/8.e3)**3)**par[4])

def BackgroundFit_f1_trigcorr_bbl(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3)))) * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f2_trigcorr_bbl(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f3_trigcorr_bbl(x, par):
	return par[0] / (1 + par[1] * (x[0] / 8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f4_trigcorr_bbl(x, par):
	if 1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2 <= 0:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f5_trigcorr_bbl(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2] * trigger_efficiency.trigger_efficiency_bbl(x[0])

def BackgroundFit_f1_trigcorr_bbh(x, par):
	return par[0] * (1. - (x[0] / 8.e3))**par[1] / ((x[0] / 8.e3)**(par[2] + par[3] * TMath.Log((x[0] / 8.e3)))) * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f2_trigcorr_bbh(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f3_trigcorr_bbh(x, par):
	return par[0] / (1 + par[1] * (x[0] / 8.e3))**par[2] * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f4_trigcorr_bbh(x, par):
	if (1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2) <= 0.:
		return 0
	elif ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) < 1.e-15:
		return 0
	else:
		return par[0] / ((1 + par[1]*x[0]/8.e3 + par[2] * (x[0]/8.e3)**2)**par[3]) * trigger_efficiency.trigger_efficiency_bbh(x[0])

def BackgroundFit_f5_trigcorr_bbh(x, par):
	return par[0] * (x[0]/8.e3)**(-1.*par[1]) * (1. - (x[0]/8.e3)**(1./3.))**par[2] * trigger_efficiency.trigger_efficiency_bbh(x[0])
