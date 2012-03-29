#!/bin/csh -f
tar zcvf tinysim.tgz cellmd eam_spe pots/ag.* pots/cu.* a.inp clusters/108.inp clusters/256.inp
tar zcvf smallsim.tgz cellmd eam_spe pots/ag.* pots/cu.* data/8k.inp
#tar zcvf largesim.tgz cellmd eam_spe pots/ag.* pots/cu.* 1m.inp
