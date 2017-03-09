import ROOT
import os
import sys
import numpy as np
import array
import math
import itertools

################# - FUNCTIONS - #################

# fill TProfiles
def FillProfiles(phi,d0_bs,tt_d0_bs):
	phi_combs	= list(itertools.combinations(phi,2))
	d0_combs	= list(itertools.combinations(d0_bs,2))
	tt_d0_combs = list(itertools.combinations(tt_d0_bs,2))
	
	ncombs = len(phi_combs)
	for i in range(ncombs):
		prof_plus_bs		.Fill(math.cos(phi_combs[i][0] + phi_combs[i][1]), d0_combs[i][0]*d0_combs[i][1] )
		prof_minus_bs		.Fill(math.cos(phi_combs[i][0] - phi_combs[i][1]), d0_combs[i][0]*d0_combs[i][1] )
		prof_plus_tt_d0_bs	.Fill(math.cos(phi_combs[i][0] + phi_combs[i][1]), tt_d0_combs[i][0]*tt_d0_combs[i][1] )
		prof_minus_tt_d0_bs	.Fill(math.cos(phi_combs[i][0] - phi_combs[i][1]), tt_d0_combs[i][0]*tt_d0_combs[i][1] )

# fit profiles and get points
def GetPoint(plus, minus):
	plus .Fit("pol1")
	minus.Fit("pol1")
	
	f_m   = minus.GetFunction("pol1")
	try: f_m   .SetLineColor(ROOT.kRed)
	except: import pdb; pdb.set_trace()
	f_m	  .SetTitle("f_minus")
	f_m	  .SetName("f_minus")
	m     = f_m.GetParameter(1)
	m_err = f_m.GetParError(1)
	
	f_p   = plus.GetFunction("pol1")
	f_p   .SetLineColor(ROOT.kBlue)
	f_p   .SetTitle("f_plus")
	f_p   .SetName("f_plus")
	p     = f_p.GetParameter(1)
	p_err = f_p.GetParError(1)
	
	try:
		x = math.sqrt(abs(m - p)) * 10000
		y = math.sqrt(abs(m + p)) * 10000
		xE = 0.5*math.sqrt( ((m_err*m_err) + (p_err * p_err))/abs(m - p) ) * 10000
		yE = 0.5*math.sqrt( ((m_err*m_err) + (p_err * p_err))/abs(m + p) ) * 10000
	except:
		import pdb; pdb.set_trace()
	
	return f_m, f_p, x, y, xE, yE

# print point to screen and file
def PrintPoint(type, x, y, xE, yE):
	
	if type:
		print           '  Type: BS'
		print >>myfile, '  Type: BS'
	else:
		print           '  Type: BS Transient Track'
		print >>myfile, '  Type: BS Transient Track'
	
	print          '    Sigma_x = %d +/- %d um' % (x,xE)
	print >>myfile,'    Sigma_x = %d +/- %d um' % (x,xE)
	print          '    Sigma_y = %d +/- %d um' % (y,yE)
	print >>myfile,'    Sigma_y = %d +/- %d um' % (y,yE)


################# - GLOBAL VARIABLES - #################
maxentries = 1000*1000*1000
limit = 0.0002

pt_cut = 1.
eta_cut = 1.
trk_hits_cut = 8
pix_hits_cut = 1

dxy_bins = 2000
dxy_min = -0.03
dxy_max = 0.03

cos_bins = 40

################# - MAIN - #################
# TEXT output file
myfile = open('comb_Fill5199_output.txt', 'w+')
print '***************** Analysing Fill5199 *****************\n'
print '------------ BEGIN OF MACRO ------------'
print >>myfile, '***************** Analysing Fill5199 *****************\n'
print >>myfile, '------------ BEGIN OF MACRO ------------'

# TChain
chain = ROOT.TChain('demo/trackTree')   # name of the tree is the argument
for i in range(1,107):
	file_path = '/gwteras/cms/store/user/fbrivio/BeamSpot/ZeroBias/crab_Fill_5199_bsWidth_good/170307_155052/0000/tracksFile_'+str(i)+'.root'
	print 'Adding: ', file_path
	print >>myfile, 'Adding: ', file_path

	chain.Add(file_path);

#chain.Add('old_files/tracksFile_10.root')
#chain.Add('tracksFile_1.root')
#chain.Add('tracksFile_106.root')
#chain.Add('tracksFile_test4.root')

print 'TChain built'
print '\t Input file: ', chain.GetName()
print >>myfile, 'TChain built'
print >>myfile, '\t Input file', chain.GetName()

# OUTPUT file
outfile = ROOT.TFile('comb_Fill5199_graphs.root','RECREATE')
print '\t Output file: ', outfile.GetName()
print >>myfile, '\t Output file: ', outfile.GetName()

# TH2D histograms
#pt_bins = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0 , 5.5, 6.0 , 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]

# TH1D per IP
h_bs 	= ROOT.TH1D("h_bs"		, "h_d0_bs"		, 100, -0.05, 0.05)
h_bs_tt	= ROOT.TH1D("h_bs_tt"	, "h_d0_bs_tt"	, 100, -0.05, 0.05)

# TProfiles
prof_plus_bs = ROOT.TProfile("prof_plus_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}+#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max)
prof_plus_bs.GetXaxis().SetTitle("cos(#phi_{1}+#phi_{2})")
prof_plus_bs.GetYaxis().SetTitle("<d_{xy}^{1} d_{xy}^{2}>")
prof_plus_bs.SetMarkerStyle(22)
prof_plus_bs.SetMarkerColor(ROOT.kBlue)
prof_minus_bs = ROOT.TProfile("prof_minus_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}-#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max)
prof_minus_bs.GetXaxis().SetTitle("cos(#phi_{1}#pm#phi_{2})")
prof_minus_bs.GetYaxis().SetTitle("<d_{xy}^{1} d_{xy}^{2}>")
prof_minus_bs.SetMarkerStyle(20)
prof_minus_bs.SetMarkerColor(ROOT.kRed)

prof_plus_tt_d0_bs = ROOT.TProfile("prof_plus_tt_d0_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}+#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max)
prof_plus_tt_d0_bs.GetXaxis().SetTitle("cos(#phi_{1}+#phi_{2})")
prof_plus_tt_d0_bs.GetYaxis().SetTitle("<d_{xy}^{1} d_{xy}^{2}>")
prof_plus_tt_d0_bs.SetMarkerStyle(22)
prof_plus_tt_d0_bs.SetMarkerColor(ROOT.kBlue)
prof_minus_tt_d0_bs = ROOT.TProfile("prof_minus_tt_d0_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}+#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max)
prof_minus_tt_d0_bs.GetXaxis().SetTitle("cos(#phi_{1}#pm#phi_{2})")
prof_minus_tt_d0_bs.GetYaxis().SetTitle("<d_{xy}^{1} d_{xy}^{2}>")
prof_minus_tt_d0_bs.SetMarkerStyle(20)
prof_minus_tt_d0_bs.SetMarkerColor(ROOT.kRed)

# TGraphs (X,Y,Xtt,Ytt)
run_graph_x = ROOT.TGraphErrors()
run_graph_x.SetName("Sigma_X")
run_graph_x.SetTitle("#sigma_{x}")
run_graph_x.GetXaxis().SetTitle("Run")
run_graph_x.GetYaxis().SetTitle("#mum")

run_graph_y = ROOT.TGraphErrors()
run_graph_y.SetName("Sigma_Y")
run_graph_y.SetTitle("#sigma_{y}")
run_graph_y.GetXaxis().SetTitle("Run")
run_graph_y.GetYaxis().SetTitle("#mum")

run_graph_x_tt = ROOT.TGraphErrors()
run_graph_x_tt.SetName("Sigma_X_tt")
run_graph_x_tt.SetTitle("#sigma_{x} - tt")
run_graph_x_tt.GetXaxis().SetTitle("Run")
run_graph_x_tt.GetYaxis().SetTitle("#mum")

run_graph_y_tt = ROOT.TGraphErrors()
run_graph_y_tt.SetName("Sigma_Y_tt")
run_graph_y_tt.SetTitle("#sigma_{y} - tt")
run_graph_y_tt.GetXaxis().SetTitle("Run")
run_graph_y_tt.GetYaxis().SetTitle("#mum")

# Check number of entries in ttree and the number of vertexes
nentries = chain.GetEntries()
#nentries = 482985501;
print '\t Number of entries  = ', nentries
print >>myfile, '\t Number of entries  = ', nentries

min_run = chain.GetMinimum("Run");
max_run = chain.GetMaximum("Run");
#min_run = 278820;
#max_run = 278822;
print '\t Run numbers: ', min_run, ' - ', max_run
print >>myfile, '\t Run numbers: ', min_run,' - ', max_run


# Loop on tracks
print ' - Begin of Track Loop - '
print >>myfile, ' - Begin of Track Loop - '

# vectors for saving tracks quantities to do combinatorial
phi		 = []
d0_bs	 = []
tt_d0_bs = []

# temporary vertexID and Run number
tmp_vtx = 1
tmp_run = min_run
point = 0

# loop on tracks
for i,evt in enumerate(chain):

	if i == maxentries: break
	
	if i%100000 == 0 :
		print "Track: ",i," - run - tmp_run: ",evt.Run," - ",tmp_run
		print >>myfile, "Track: ",i," - run - tmp_run: ",evt.Run," - ",tmp_run
	
	# check when the run changes or this is the last track => new point in graph
	if evt.Run == tmp_run and i != nentries-1:

		# Track selection (pt>1 - |eta|<1 - 8Hits - 1PixHit)
		if evt.pt < pt_cut or abs(evt.eta) > eta_cut or evt.Track_HITs < trk_hits_cut or evt.Pix_HITs < pix_hits_cut : continue
			
		#print '\tTrack:', i

		# if from same vertex => save interesting quantities
		if evt.VtxID == tmp_vtx:
			phi.append(evt.phi)
			d0_bs.extend([evt.d0_bs])
			tt_d0_bs.extend([evt.tt_d0_bs])
		
		# if from different vertex:
		#	- fill the TProfiles
		#	- clear the tmp_arrays
		#	- insert new value in tmp_arrays and in tmp_vtx
		else:
			if len(phi) >1:
				FillProfiles(phi,d0_bs,tt_d0_bs)
			
			phi		 = []
			d0_bs	 = []
			tt_d0_bs = []

			phi.append(evt.phi)
			d0_bs.extend([evt.d0_bs])
			tt_d0_bs.extend([evt.tt_d0_bs])

			tmp_vtx = evt.VtxID

	# new run:
	#  - fit tprofiles
	#  - set point in tgraphs
	#  - save objects in root file
	#  - clear all
	else:
	
		# get point coordinates and fit functions
		f_m, f_p, x, y, xE, yE					 = GetPoint(prof_plus_bs		, prof_minus_bs)
		f_m_tt, f_p_tt, x_tt, y_tt, xE_tt, yE_tt = GetPoint(prof_plus_tt_d0_bs	, prof_minus_tt_d0_bs)
		
		# print point
		print           ' --- NEW POINT, Run: %d ---' % (tmp_run)
		print >>myfile, ' --- NEW POINT, Run: %d ---' % (tmp_run)
		PrintPoint(1, x, y, xE, yE)				# to screen
		PrintPoint(0, x_tt, y_tt, xE_tt, yE_tt)	# to file
		print           ' --------------------------'
		print >>myfile, ' --------------------------'
		
		# set points in tgraph
		run_graph_x.SetPoint		(point, tmp_run, x)
		run_graph_x.SetPointError	(point, 0.     , xE)
		run_graph_y.SetPoint		(point, tmp_run, y)
		run_graph_y.SetPointError	(point, 0.     , yE)

		run_graph_x_tt.SetPoint		(point, tmp_run, x_tt)
		run_graph_x_tt.SetPointError(point, 0.     , xE_tt)
		run_graph_y_tt.SetPoint		(point, tmp_run, y_tt)
		run_graph_y_tt.SetPointError(point, 0.     , yE_tt)
		
		# save to file
		outfile.cd()
		outfile.mkdir(str(tmp_run))
		outfile.cd(str(tmp_run))
		prof_plus_bs.Write()
		prof_minus_bs.Write()
		prof_plus_tt_d0_bs.Write()
		prof_minus_tt_d0_bs.Write()
		f_m.Write()
		f_p.Write()
		f_m_tt.Write()
		f_p_tt.Write()
		outfile.cd()
		
		# clear all
		prof_minus_bs		.Reset()
		prof_plus_bs		.Reset()
		prof_plus_tt_d0_bs	.Reset()
		prof_minus_tt_d0_bs	.Reset()

		phi		 = []
		d0_bs	 = []
		tt_d0_bs = []
		
		phi.append(evt.phi)
		d0_bs.extend([evt.d0_bs])
		tt_d0_bs.extend([evt.tt_d0_bs])
		
		tmp_vtx = evt.VtxID
		tmp_run = evt.Run
		point = point+1
		
#import pdb; pdb.set_trace()

outfile.cd()
run_graph_x.Write()
run_graph_y.Write()
run_graph_x_tt.Write()
run_graph_y_tt.Write()
outfile.Close();
print "------------ END OF MACRO ------------"
print >>myfile,"------------ END OF MACRO ------------"
myfile.close();







'''
for i in range(chain.GetEntries()):
	
	chain.GetEntry(i)
	
	if i == maxentries : break
	
	if i%100000 == 0 :
		print "Track: ",i," - run - tmp_run: ",Runi," - ",tmp_run
		print >>myfile, "Track: ",i," - run - tmp_run: ",Runi," - ",tmp_run

	errIP_pt.Fill(ipt[0],itt_d0_err_bs[0])
	errIP_eta.Fill(ieta[0],itt_d0_err_bs[0])
	errIP_phi.Fill(iphi[0],itt_d0_err_bs[0])


	Runi  = iRun[0]

	# Track selection (pt>1 - |eta|<1 - 8Hits - 1PixHit)
	if ipt[0] < pt_cut or abs(ieta[0]) > eta_cut or iTrack_HITs[0] < trk_hits_cut or iPix_HITs[0] < pix_hits_cut : continue

	after_cuts += 1

	# Temporary variables
	tmp_vtx = iVtxID[0]
	tmp_phi = iphi[0]
	tmp_d0_bs = id0_bs[0]
	tmp_tt_d0_bs = itt_d0_bs[0]

	h_bs.Fill(id0_bs[0])

	# Second loop on tracks
	for j in range (i+1,1000):
	
		chain.GetEntry(j)
	
		# check if the vetex is the same
		if iVtxID[0] > tmp_vtx: break
		
		# Track selection (pt>1 - |eta|<1 - 8Hits - 1PixHit)
		if ipt[0] < pt_cut or abs(ieta[0]) > eta_cut or iTrack_HITs[0] < trk_hits_cut or iPix_HITs[0] < pix_hits_cut : continue

		# Fill TProfiles
		prof_plus_bs.Fill		(math.cos(iphi[0] + tmp_phi), id0_bs[0]*tmp_d0_bs 	);
		prof_minus_bs.Fill		(math.cos(iphi[0] - tmp_phi), id0_bs[0]*tmp_d0_bs 	);
		prof_plus_tt_d0_bs.Fill	(math.cos(iphi[0] + tmp_phi), itt_d0_bs[0]*tmp_tt_d0_bs);
		prof_minus_tt_d0_bs.Fill(math.cos(iphi[0] - tmp_phi), itt_d0_bs[0]*tmp_tt_d0_bs);
		
		if iRun[0] == min_run: first_run += 1


	if tmp_run == Runi:
		tmp_vtx      = 0
		tmp_phi      = 0.
		tmp_d0_bs    = 0.
		tmp_tt_d0_bs = 0.
		continue
	else:
		print '\n - NEW POINT IN GRAPH, run - tmp_run: ',Runi,' - ',tmp_run
		print >>myfile,'\n - NEW POINT IN GRAPH, run - tmp_run: ',Runi,' - ',tmp_run
		try:
			prof_minus_bs.SetMinimum(-1.*limit)
			prof_minus_bs.SetMaximum(limit)
			prof_plus_bs .SetMinimum(-1.*limit)
			prof_plus_bs .SetMaximum(limit)
			prof_minus_bs.Fit("pol1")
			prof_plus_bs .Fit("pol1")
			
			f_m_bs = prof_minus_bs.GetFunction("pol1");
			f_p_bs = prof_plus_bs.GetFunction("pol1");
			f_p_bs.SetLineColor(ROOT.kBlue);
			f_p_bs.SetLineColor(ROOT.kRed);
			m_bs = f_m_bs.GetParameter(1)
			m_bs_err = f_m_bs.GetParError(1);
			p_bs = f_p_bs.GetParameter(1)
			p_bs_err = f_p_bs.GetParError(1);
			sigma_x_err = 0.5*math.sqrt(abs( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs - p_bs) ));
			sigma_y_err = 0.5*math.sqrt(abs( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs + p_bs) ));


			print "Sigma_x = ",math.sqrt(m_bs - p_bs) * 10000," +/- ",sigma_x_err * 10000," um"
			print "Sigma_y = ",math.sqrt(m_bs + p_bs) * 10000," +/- ",sigma_y_err * 10000," um"
			print >>myfile, "Sigma_x = ",math.sqrt(m_bs - p_bs) * 10000," +/- ",sigma_x_err * 10000," um"
			print >>myfile, "Sigma_y = ",math.sqrt(m_bs + p_bs) * 10000," +/- ",sigma_y_err * 10000," um"


			print "\t Run:",tmp_run
			print "\t x:"<<math.sqrt(m_bs - p_bs)*10000
			print "\t y:"<<math.sqrt(m_bs + p_bs)*10000
			print >>myfile, "\t Run:",tmp_run
			print >>myfile, "\t x:"<<math.sqrt(m_bs - p_bs)*10000
			print >>myfile, "\t y:"<<math.sqrt(m_bs + p_bs)*10000
			run_graph_x.SetPoint(point, tmp_run, math.sqrt(m_bs - p_bs)*10000);
			run_graph_x.SetPointError(point,0., sigma_x_err*10000);
			run_graph_y.SetPoint(point, tmp_run, math.sqrt(m_bs + p_bs)*10000);
			run_graph_y.SetPointError(point,0., sigma_y_err*10000);
			
			prof_minus_bs.Reset();
			prof_plus_bs.Reset();

			tmp_run = Runi;
			point = poin+1;
		
			#Clear temporary variables
			tmp_vtx      = 0;
			tmp_phi      = 0.;
			tmp_d0_bs    = 0.;
			tmp_tt_d0_bs = 0.;
		except:
			print "Unexpected error:", sys.exc_info()[0]
			raise

print " - End of Track Loop - \n";
print "\t Number of events after cuts = ",after_cuts
print "\t \t First Run events : ",first_run
print "\t \t Last Run events  : ",last_run
print >>myfile, " - End of Track Loop - \n";
print >>myfile, "\t Number of events after cuts = ",after_cuts
print >>myfile, "\t \t First Run events : ",first_run
print >>myfile, "\t \t Last Run events  : ",last_run

# Fit the TProfiles
prof_minus_bs.SetMinimum(-1.*limit)
prof_minus_bs.SetMaximum(limit)
prof_plus_bs .SetMinimum(-1.*limit)
prof_plus_bs .SetMaximum(limit)

prof_minus_tt_d0_bs.SetMinimum(-1.*limit)
prof_minus_tt_d0_bs.SetMaximum(limit)
prof_plus_tt_d0_bs .SetMinimum(-1.*limit)
prof_plus_tt_d0_bs .SetMaximum(limit)

prof_minus_bs.Fit("pol1")
prof_plus_bs .Fit("pol1");
prof_minus_tt_d0_bs.Fit("pol1");
prof_plus_tt_d0_bs .Fit("pol1",'','',-1.,0.1);

# Retrieve fit parameters
f_m_bs = prof_minus_bs.GetFunction("pol1");
f_p_bs = prof_plus_bs.GetFunction("pol1");
f_p_bs.SetLineColor(ROOT.kBlue);
f_p_bs.SetLineColor(ROOT.kRed);
m_bs = f_m_bs.GetParameter(1)
m_bs_err = f_m_bs.GetParError(1)
m_chi2 = (f_m_bs.GetChisquare())/(f_m_bs.GetNDF());
p_bs = f_p_bs.GetParameter(1)
p_bs_err = f_p_bs.GetParError(1);
p_chi2 = (f_p_bs.GetChisquare())/(f_p_bs.GetNDF());
sigma_x_err = 0.5#*math.sqrt( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs - p_bs) );
sigma_y_err = 0.5#*math.sqrt( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs + p_bs) );

import pdb; pdb.set_trace()

print "-------- LAST ONE ---------\n";
print "\n - NEW POINT IN GRAPH, run - tmp_run: ",Runi," - ",tmp_run
print "Sigma_x = ",math.sqrt(m_bs - p_bs) * 10000," +/- ",sigma_x_err * 10000," um"
print "Sigma_y = ",math.sqrt(m_bs + p_bs) * 10000," +/- ",sigma_y_err * 10000," um"
print >>myfile, "-------- LAST ONE ---------\n";
print >>myfile, "\n - NEW POINT IN GRAPH, run - tmp_run: ",Runi," - ",tmp_run
print >>myfile, "Sigma_x = ",math.sqrt(m_bs - p_bs) * 10000," +/- ",sigma_x_err * 10000," um"
print >>myfile, "Sigma_y = ",math.sqrt(m_bs + p_bs) * 10000," +/- ",sigma_y_err * 10000," um"
run_graph_x.SetPoint(point, tmp_run, math.sqrt(m_bs - p_bs)*10000);
run_graph_x.SetPointError(point,0., sigma_x_err*10000);
run_graph_y.SetPoint(point, tmp_run, math.sqrt(m_bs + p_bs)*10000);
run_graph_y.SetPointError(point,0., sigma_y_err*10000);

f_m_tt = prof_minus_tt_d0_bs.GetFunction("pol1");
f_p_tt = prof_plus_tt_d0_bs .GetFunction("pol1");
m_tt = f_m_tt.GetParameter(1);
p_tt = f_p_tt.GetParameter(1);


print "\n*************** - FIT RESULTS - ***************"
print "	- bs - "
print "m_bs = ",m_bs," +/- ",m_bs_err
print "m_chi2 = ",f_m_bs.GetChisquare(),"/",f_m_bs.GetNDF()," = ",m_chi2
print "p_bs = ",p_bs," +/- ",p_bs_err
print "p_chi2 = ",f_p_bs.GetChisquare(),"/",f_p_bs.GetNDF()," = ",p_chi2
print "Sigma_x^2 = ",m_bs - p_bs
print "Sigma_y^2 = ",m_bs + p_bs
print "Sigma_x = ",math.sqrt(m_bs - p_bs) * 10000," +/- ",sigma_x_err * 10000," um"
print "Sigma_y = ",math.sqrt(m_bs + p_bs) * 10000," +/- ",sigma_y_err * 10000," um"
print "\n";
print "	- tt_bs - "
print "Sigma_x^2 = ",m_tt - p_tt
print "Sigma_y^2 = ",m_tt + p_tt
print "Sigma_x = ",math.sqrt(m_tt - p_tt) * 10000," um"
print "Sigma_y = ",math.sqrt(m_tt + p_tt) * 10000," um"
print "*************************************************"

print >>myfile, "\n*************** - FIT RESULTS - ***************"
print >>myfile, "	- bs - "
print >>myfile, "m_bs = ",m_bs," +/- ",m_bs_err
print >>myfile, "m_chi2 = ",f_m_bs.GetChisquare(),"/",f_m_bs.GetNDF()," = ",m_chi2
print >>myfile, "p_bs = ",p_bs," +/- ",p_bs_err
print >>myfile, "p_chi2 = ",f_p_bs.GetChisquare(),"/",f_p_bs.GetNDF()," = ",p_chi2
print >>myfile, "Sigma_x^2 = ",m_bs - p_bs
print >>myfile, "Sigma_y^2 = ",m_bs + p_bs
print >>myfile, "Sigma_x = ",math.sqrt(m_bs - p_bs) * 10000," +/- ",sigma_x_err * 10000," um"
print >>myfile, "Sigma_y = ",math.sqrt(m_bs + p_bs) * 10000," +/- ",sigma_y_err * 10000," um"
print >>myfile, "\n";
print >>myfile, "	- tt_bs - "
print >>myfile, "Sigma_x^2 = ",m_tt - p_tt
print >>myfile, "Sigma_y^2 = ",m_tt + p_tt
print >>myfile, "Sigma_x = ",math.sqrt(m_tt - p_tt) * 10000," um"
print >>myfile, "Sigma_y = ",math.sqrt(m_tt + p_tt) * 10000," um"
print >>myfile, "*************************************************"

# Projection of TH2D of errIP
proj_errIP_pt = errIP_pt.ProjectionX();
proj_errIP_phi = errIP_phi.ProjectionX();
proj_errIP_eta = errIP_eta.ProjectionX();


# Close file and save the results
outfile.cd();

proj_errIP_pt	. Write();
proj_errIP_eta	. Write();
proj_errIP_phi	. Write();

h_bs			. Write();

prof_plus_bs		. Write();
prof_minus_bs		. Write();
prof_plus_tt_d0_bs 	. Write();
prof_minus_tt_d0_bs	. Write();

run_graph_x		. Write();
run_graph_y		. Write();

outfile			. Close();

print "------------ END OF MACRO ------------"
print >>myfile,"------------ END OF MACRO ------------"
myfile.close();


'''











