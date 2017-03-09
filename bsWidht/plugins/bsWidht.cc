// -*- C++ -*-
//
// Package:    bsWidth/bsWidht
// Class:      bsWidht
// 
/**\class bsWidht bsWidht.cc bsWidth/bsWidht/plugins/bsWidht.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco Brivio
//         Created:  Fri, 17 Feb 2017 20:26:05 GMT
//
//


// system include files
#include <memory>
#include <algorithm>
#include <vector>
#include "TRandom.h"
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class bsWidht : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit bsWidht(const edm::ParameterSet&);
      ~bsWidht();

	  void emptyAll();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	  static bool mysorter (reco::Track i, reco::Track j) { return (i.pt () > j.pt()); }

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
	
	  edm::Service<TFileService> outfile_;
	  TTree* tree_;
	
	  edm::InputTag pvsTag_;
	  edm::EDGetTokenT<reco::VertexCollection> pvsToken_;
	
	  edm::InputTag tracksTag_;
	  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
	
	  edm::InputTag bsTag_;
	  edm::EDGetTokenT<reco::BeamSpot> bsToken_;
	
	  // Variables for branches
	  int Run_;
	  int Lumi_;
	  int Event_;
	  int Collision_ = 0;
	  int VtxID_     = 0;
	  int vtx_ntrks_;
	  double x_PV_;
	  double y_PV_;
	  double z_PV_;

	  double pt_;
	  double eta_;
	  double phi_;
	  double chi2_;
	  double IP_;
	  double d0_;
	  double d0_bs_;
	  double d0_xyz_;
	  double dZ_bs_;
	  double tt_d0_bs_;
	  double tt_d0_err_bs_;
	  double tt_z0_bs_;
	  int Pix_HITs_;
	  int Strip_HITs_;
	  int Track_HITs_;
	  bool high_quality_;
	
	  // beamSpot object
	  reco::BeamSpot beamSpot;
	



      // ----------member data ---------------------------
};

//
// constructors and destructor
//
bsWidht::bsWidht(const edm::ParameterSet& iConfig):
pvsTag_      (iConfig.getParameter<edm::InputTag>("vtxCollection")),
pvsToken_    (consumes<reco::VertexCollection>(pvsTag_)),
tracksTag_   (iConfig.getParameter<edm::InputTag>("trackCollection")),
tracksToken_ (consumes<reco::TrackCollection>(tracksTag_)),
bsTag_		 (iConfig.getParameter<edm::InputTag>("beamSpot")),
bsToken_     (consumes<reco::BeamSpot>(bsTag_))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
}

bsWidht::~bsWidht()
{
}


//
// member functions
//

// ------------ method called for each event  ------------
void
bsWidht::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
	
   emptyAll();
	
   //Fill general info
	Run_		= iEvent.id().run();
	Lumi_		= iEvent.id().luminosityBlock();
	Event_		= iEvent.id().event();
	Collision_  = Collision_ +1;
	
	// Handles
	edm::ESHandle<TransientTrackBuilder> theB;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);
	
	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(pvsToken_, vertices);
	const reco::VertexCollection pvtx  = *(vertices.product())  ;
	
	edm::Handle<reco::TrackCollection> tracks;
	iEvent.getByToken(tracksToken_, tracks);
	
	edm::Handle<reco::BeamSpot> bs;
	iEvent.getByToken(bsToken_, bs);
	
	for (reco::VertexCollection::const_iterator pvIt = pvtx.begin(); pvIt!=pvtx.end(); pvIt++)
	{
		reco::Vertex iPV = *pvIt;
		if (iPV.isFake()) continue;
		VtxID_ += 1;
		//std::cout<<"Vertex n: "<<VtxID_<<std::endl;
		x_PV_ 	 	= iPV.x();
		y_PV_ 	 	= iPV.y();
		z_PV_		= iPV.z();
		
		// get all valid tracks
		reco::TrackCollection allTracks;
		reco::Vertex::trackRef_iterator trki;
		for (trki  = iPV.tracks_begin(); trki != iPV.tracks_end(); ++trki)
		{
			if (trki->isNonnull())
			{
				reco::TrackRef trk_now(tracks, (*trki).key());
				allTracks.push_back(*trk_now);
			}
		}
		
		// order with decreasing pt
		//std::sort (allTracks.begin(), allTracks.end(), mysorter);
		
		// count good tracks per vertex
		int good_tracks = 0;
		uint ntrks = allTracks.size();
		for (uint tracksIt = 0 ;  tracksIt < ntrks; tracksIt = tracksIt+1)
		{
			reco::Track firstTrk = allTracks.at(tracksIt);
			if (firstTrk.numberOfValidHits()>=8 && firstTrk.pt()>=1. && firstTrk.normalizedChi2()<10.
				&& firstTrk.normalizedChi2()>0. && firstTrk.quality(reco::TrackBase::highPurity)
				&& (std::abs(firstTrk.eta()) < 1.) && firstTrk.hitPattern().numberOfValidPixelHits() >=1)
			{
				// increase number of good tracks per vertex
				good_tracks += 1;
			}
		}
		vtx_ntrks_ = good_tracks;
		//std::cout << "\t good tracks:" << vtx_ntrks_ << std::endl;
		
		// loop on tracks to save the quantities
		for (uint tracksIt = 0 ;  tracksIt < ntrks; tracksIt = tracksIt+1)
		{
			reco::Track firstTrk = allTracks.at(tracksIt);
			if (firstTrk.numberOfValidHits()>=8 && firstTrk.pt()>=1. && firstTrk.pt()<5000. && firstTrk.normalizedChi2()<10.
				&& firstTrk.normalizedChi2()>0. && firstTrk.quality(reco::TrackBase::highPurity)
				&& (std::abs(firstTrk.eta()) < 1.) && firstTrk.hitPattern().numberOfValidPixelHits() >=1)
			{
				// track varaibles
				pt_				= firstTrk.pt();
				eta_			= firstTrk.eta();
				chi2_			= firstTrk.normalizedChi2();
				phi_			= firstTrk.phi();
				IP_				= firstTrk.dxy();
				Pix_HITs_		= firstTrk.hitPattern().numberOfValidPixelHits();
				Strip_HITs_		= firstTrk.hitPattern().numberOfValidStripHits();
				Track_HITs_		= firstTrk.numberOfValidHits();
				high_quality_	= firstTrk.quality(reco::TrackBase::highPurity);
				
				// impact parameters
				math::XYZPoint zero_point(0.,0.,0.);
				d0_ = -1.* firstTrk.dxy(zero_point);
				
				math::XYZPoint PV_point(x_PV_, y_PV_, z_PV_);
				d0_xyz_ = -1.* firstTrk.dxy(PV_point);
				
				beamSpot = *bs;
				math::XYZPoint BSpoint(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
				d0_bs_ = -1.* firstTrk.dxy(BSpoint);
				dZ_bs_ = -1 * firstTrk.dz(BSpoint);
				
				const reco::TransientTrack trans_track = (*theB).build(firstTrk);
				GlobalPoint BSvert(beamSpot.x0(),beamSpot.y0(), beamSpot.z0());
				TrajectoryStateClosestToPoint  traj = trans_track.trajectoryStateClosestToPoint(BSvert);
				tt_d0_bs_ = traj.perigeeParameters().transverseImpactParameter();
				tt_d0_err_bs_ = traj.perigeeError().transverseImpactParameterError();
				tt_z0_bs_ = traj.perigeeParameters().longitudinalImpactParameter();
			}
			else
			{
				pt_				= -99.;
				eta_			= -99.;
				chi2_			= -99.;
				phi_			= -99.;
				IP_				= -99.;
				Pix_HITs_		= -99.;
				Strip_HITs_		= -99.;
				Track_HITs_		= -99.;
				high_quality_	= -99.;
				d0_				= -99.;
				d0_xyz_			= -99.;
				d0_bs_			= -99.;
				dZ_bs_			= -99.;
				tt_d0_bs_		= -99.;
				tt_d0_err_bs_	= -99.;
				tt_z0_bs_		= -99.;
			}
			
			//Fill tree
			tree_ -> Fill();
			
		} // Track loop to fill varaibles and Impact parameters

		
	} // Vertex Loop

} // end function Analyze

// ------------ method called once each job just before starting event loop  ------------
void 
bsWidht::beginJob()
{
	// Create tree
	tree_ = outfile_-> make<TTree>("trackTree","trackTree");
	
	// TTree branches
	tree_->Branch("Run"       	,&Run_		);
	tree_->Branch("Lumi"		,&Lumi_		);
	tree_->Branch("Event"		,&Event_	);
	tree_->Branch("Collision"	,&Collision_);
	
	tree_->Branch("VtxID"		,&VtxID_	);
	tree_->Branch("x_PV"	   	,&x_PV_		);
	tree_->Branch("y_PV"	   	,&y_PV_		);
	tree_->Branch("z_PV"	   	,&z_PV_		);
	tree_->Branch("vtx_ntrks"	,&vtx_ntrks_);
	
	tree_->Branch("pt"			,&pt_			);
	tree_->Branch("eta"			,&eta_			);
	tree_->Branch("phi"			,&phi_			);
	tree_->Branch("chi2"		,&chi2_			);
	tree_->Branch("dxy"       	,&IP_			);
	tree_->Branch("d0"			,&d0_			);
	tree_->Branch("d0_bs"		,&d0_bs_		);
	tree_->Branch("d0_xyz"		,&d0_xyz_		);
	tree_->Branch("dZ_bs"		,&dZ_bs_		);
	tree_->Branch("tt_d0_bs"	,&tt_d0_bs_		);
	tree_->Branch("tt_d0_err_bs",&tt_d0_err_bs_	);
	tree_->Branch("tt_z0_bs"	,&tt_z0_bs_		);
	
	tree_->Branch("Pix_HITs"	,&Pix_HITs_		);
	tree_->Branch("Strip_HITs"	,&Pix_HITs_		);
	tree_->Branch("Track_HITs"	,&Track_HITs_	);
	tree_->Branch("high_quality",&high_quality_	);
	


}

// ------------ method called once each job just after ending the event loop  ------------
void 
bsWidht::endJob() 
{
	tree_ = tree_->CopyTree("eta>-10. && d0_bs<100. && d0_bs>-100. && tt_d0_bs<100. && tt_d0_bs>-100.");
	gDirectory->Purge();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
bsWidht::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// ----------- empy all variables ----------
void bsWidht::emptyAll()
{
	Run_ = 0;
	Lumi_= 0;
	Event_= 0;
	VtxID_ = 0;
	vtx_ntrks_ = 0;
	x_PV_ = 0.;
	y_PV_ = 0.;
	z_PV_ = 0.;
	pt_ = 0.;
	eta_ = 0.;
	phi_ = 0.;
	chi2_ = 0.;
	IP_ = 0.;
	d0_ = 0.;
	d0_bs_ = 0.;
	d0_xyz_ = 0.;
	dZ_bs_ = 0.;
	tt_d0_bs_ = 0.;
	tt_d0_err_bs_ = 0.;
	tt_z0_bs_ = 0.;
	Pix_HITs_ = 0;
	Strip_HITs_ = 0;
	Track_HITs_ = 0;
	high_quality_ = false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(bsWidht);
