import os
import sys
import pickle
import math

cL_uptype = 0.3442*2
cR_uptype = -0.1558*2
cL_downtype = -0.4221*2
cR_downtype = 0.0779*2
cuuL = cL_uptype
cuuR = cR_uptype
cccL = cL_uptype
cccR = cR_uptype
cttL = cL_uptype
cttR = cR_uptype
cddL = cL_downtype
cddR = cR_downtype
cssL = cL_downtype
cssR = cR_downtype
cbbL = cL_downtype
cbbR = cR_downtype
aEW = 1. / 127.9 
ee = 2*math.sqrt(aEW)*math.sqrt(math.pi)
MB = 4.7
MT = 172.44
MW = 80.
MZ = 91.2
sw2 = 1 - MW**2/MZ**2
sw = math.sqrt(sw2)
cw = math.sqrt(1. - sw2)
def get_partial_width(decay, MZP=1000.):
	if decay == "bb":
		return (((-6*cbbL**2*ee**2*MB**2)/(cw**2*sw**2) + (36*cbbL*cbbR*ee**2*MB**2)/(cw**2*sw**2) - (6*cbbR**2*ee**2*MB**2)/(cw**2*sw**2) + (6*cbbL**2*ee**2*MZP**2)/(cw**2*sw**2) + (6*cbbR**2*ee**2*MZP**2)/(cw**2*sw**2))*math.sqrt(-4*MB**2*MZP**2 + MZP**4))/(48.*math.pi*abs(MZP)**3)
	elif decay == "cc":
		return (MZP**2*((6*cccL**2*ee**2*MZP**2)/(cw**2*sw**2) + (6*cccR**2*ee**2*MZP**2)/(cw**2*sw**2)))/(48.*math.pi*abs(MZP)**3)
	elif decay == "dd":
		return (MZP**2*((6*cddL**2*ee**2*MZP**2)/(cw**2*sw**2) + (6*cddR**2*ee**2*MZP**2)/(cw**2*sw**2)))/(48.*math.pi*abs(MZP)**3)
	elif decay == "ss":
		return (MZP**2*((6*cssL**2*ee**2*MZP**2)/(cw**2*sw**2) + (6*cssR**2*ee**2*MZP**2)/(cw**2*sw**2)))/(48.*math.pi*abs(MZP)**3)
	elif decay == "tt":
		if MZP > 2 * MT:
			return (((-6*cttL**2*ee**2*MT**2)/(cw**2*sw**2) + (36*cttL*cttR*ee**2*MT**2)/(cw**2*sw**2) - (6*cttR**2*ee**2*MT**2)/(cw**2*sw**2) + (6*cttL**2*ee**2*MZP**2)/(cw**2*sw**2) + (6*cttR**2*ee**2*MZP**2)/(cw**2*sw**2))*math.sqrt(-4*MT**2*MZP**2 + MZP**4))/(48.*math.pi*abs(MZP)**3)
		else:
			return 0.
	elif decay == "uu":
		return (MZP**2*((6*cuuL**2*ee**2*MZP**2)/(cw**2*sw**2) + (6*cuuR**2*ee**2*MZP**2)/(cw**2*sw**2)))/(48.*math.pi*abs(MZP)**3)
	else:
		print "[get_partial_width] ERROR : decay must be in [uu, dd, cc, ss, tt, bb]"
		sys.exit(1)

def get_total_width(MZP=1000.):
	return sum([get_partial_width(x, MZP) for x in ["uu", "dd", "cc", "ss", "tt", "bb"]])

#masses = [325]
#masses.extend(range(350, 1250, 50))
masses = [400]
basedir = "/home/dryu/Dijets/data/EightTeeEeVeeBee/ZPrime/MG5_aMC_v2_5_5/Zp2b"
final_states = ["u u~", "d d~", "c c~", "s s~", "b b~"]
#final_states = ["b b~"]
mt = 172.44
zeta_jj = {}
zeta_bb = {}
sigmatilde = {}
xs = {}
for mass in masses:
	zeta_jj[mass] = {}
	zeta_bb[mass] = {}
	sigmatilde[mass] = {}
	xs[mass] = {}
	for name, initial_states in {"u":["u u~"], "d":["d d~"], "pp":["u u~", "d d~"]}.iteritems():
		tag = str(mass) + "_" + name
		width_file = open("{}/{}_{}/Cards/param_card.dat".format(basedir, mass, name), 'r')
		#width = 0
		#for line in width_file:
		#	if "DECAY 32" in line:
		#		width = float(line.split()[2])
		#		break
		#if width == 0:
		#	print "ERROR : Couldn't find width."
		#	sys.exit(1)
		width = get_total_width(mass)
		width_file.close()
		xs_file = open("{}/{}_{}/cross_section.txt".format(basedir, mass, name), 'r')
		xs[mass][name] = 0
		for line in xs_file:
			if "run_01" in line:
				xtra_fac = 1.
				if name in ["u", "d"]:
					xtra_fac = 2.
				xs[mass][name] = float(line.split()[2]) * xtra_fac # Factor of two from MG explicitly taking u from PDF1 and u~ from PDF2
		if xs[mass][name] == 0:
			print "ERROR : Couldn't find xs"
			sys.exit(1)

		if name == "pp":
			BR_initial = (get_partial_width("uu", mass) + get_partial_width("dd", mass)) / get_total_width(mass)
		elif name == "u":
			BR_initial = get_partial_width("uu", mass) / get_total_width(mass)
		elif name == "d":
			BR_initial = get_partial_width("dd", mass) / get_total_width(mass)
		BR_final_jj = 1. # All quarks
		zeta_jj[mass][name] = width / mass * BR_initial * BR_final_jj
		BR_final_bb = get_partial_width("bb", mass) / get_total_width(mass)
		zeta_bb[mass][name] = width / mass * BR_initial * BR_final_bb
		sigmatilde[mass][name] = xs[mass][name] / zeta_jj[mass][name]
		print "sigmatilde[{}][{}] = {} / {} = {}".format(mass, name, xs[mass][name], zeta_jj[mass][name], sigmatilde[mass][name])
		print "\twidth = {}".format(width)
pickle.dump(zeta_jj, open("{}/zeta2_jj.pkl".format(basedir), "wb"))
pickle.dump(zeta_bb, open("{}/zeta2_bb.pkl".format(basedir), "wb"))
pickle.dump(sigmatilde, open("{}/sigmatilde2.pkl".format(basedir), "wb"))
pickle.dump(xs, open("{}/xs2.pkl".format(basedir), "wb"))


