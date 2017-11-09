import os
import sys
import pickle
import math

masses = [325]
masses.extend(range(350, 1250, 50))
masses = [400]
basedir = "/home/dryu/Dijets/data/EightTeeEeVeeBee/ZPrime/MG5_aMC_v2_5_5/Zp1b"
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
		#width_file = open("{}/{}_{}/Cards/param_card.dat".format(basedir, mass, name), 'r')
		#width = 0
		#for line in width_file:
		#	if "DECAY 10030" in line:
		#		width = float(line.split()[2])
		#		break
		#if width == 0:
		#	print "ERROR : Couldn't find width in file {}.".format("{}/{}_{}/Cards/param_card.dat".format(basedir, mass, name))
		#	sys.exit(1)
		#width_file.close()
		width = (0.75**2 * mass * (5. + math.sqrt(max(0., 1 - (4*mt**2)/mass**2))*(1 + 2*mt**2/mass**2)))/(144.*math.pi)

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

		if mass > 2 * mt:
			tpart = (1. - 4*mt**2/mass**2)**0.5 * (1 + 2*mt**2/mass**2)
		else:
			tpart = 0.
		BR_initial = len(initial_states) / (5 + tpart)
		BR_final_jj = 1. # All quarks
		zeta_jj[mass][name] = width / mass * BR_initial * BR_final_jj
		BR_final_bb = 1. / (5. + tpart)
		zeta_bb[mass][name] = width / mass * BR_initial * BR_final_bb
		sigmatilde[mass][name] = xs[mass][name] / zeta_jj[mass][name]
		print xs[mass][name]
		print "sigmatilde[{}][{}] = {} / {} = {}".format(mass, name, xs[mass][name], zeta_jj[mass][name], sigmatilde[mass][name])
		print "\twidth = {}".format(width)
pickle.dump(zeta_jj, open("{}/zeta_jj.pkl".format(basedir), "wb"))
pickle.dump(zeta_bb, open("{}/zeta_bb.pkl".format(basedir), "wb"))
pickle.dump(sigmatilde, open("{}/sigmatilde.pkl".format(basedir), "wb"))
pickle.dump(xs, open("{}/xs.pkl".format(basedir), "wb"))


