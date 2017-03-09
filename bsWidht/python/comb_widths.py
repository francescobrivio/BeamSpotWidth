import ROOT
import os
import sys
import numpy as np
import array
import math

maxentries = 10000000000000000

# Global variables
pt_cut = 1.
eta_cut = 1.
trk_hits_cut = 8
pix_hits_cut = 1

dxy_bins = 2000
dxy_min = -0.03
dxy_max = 0.03

cos_bins = 40


# TChain
chain = ROOT.TChain('demo/trackTree')   # name of the tree is the argument
chain.Add('tracksFile_test4.root')


# vectors for saving tracks quantities to do combinatorial
phi		= []
d0_bs		= []
tt_d0_bs	= []

# temporary vertex id variable
tmp_vtx = 1

# loop on tracks
for i,evt in enumerate(chain):

	# Track selection (pt>1 - |eta|<1 - 8Hits - 1PixHit)
	if evt.pt < pt_cut or abs(evt.eta) > eta_cut or evt.Track_HITs < trk_hits_cut or evt.Pix_HITs < pix_hits_cut : continue

	print 'Track:', i

	# save interesting quantities
	if (evt.VtxID == tmp_vtx):
		phi.append([evt.phi])
		d0_bs.extend([evt.d0_bs])
		tt_d0_bs.extend([evt.tt_d0_bs])
	else:
		import pdb; pdb.set_trace()















