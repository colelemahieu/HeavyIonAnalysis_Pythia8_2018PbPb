// -*- C++ -*-
//
// Package:    ZDCTreeProducer
// Class:      ZDCTreeProducer
//
/**\class ZDCTreeProducer ZDCTreeProducer.cc CmsHi/ZDCTreeProducer/src/ZDCTreeProducer.cc
   Description: [one line class summary]
   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Yetkin Yilmaz
// Modified: Frank Ma, Yen-Jie Lee
//         Created:  Tue Sep  7 11:38:19 EDT 2010
// $Id: RecHitTreeProducer.cc,v 1.27 2013/01/22 16:36:27 yilmaz Exp $
//
//

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibCalorimetry/HcalAlgos/interface/HcalPulseShapes.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/METReco/interface/HcalCaloFlagLabels.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "DataFormats/HcalDigi/interface/HcalQIESample.h"
//#include "QWAna/QWZDC2018RecHit/src/QWZDC2018Helper.h"
#include "HeavyIonsAnalysis/ZDCAnalysis/src/QWZDC2018Helper.h"


#include "TTree.h"
#include "TNtuple.h"

#define MAXHITS 100000

struct MyZDCRecHit{
  int n;
  float  e[18];
  int    zside[18];
  int    section [18];
  int    channel[18];
  int    saturation[18];
};

struct MyZDCDigi{
  int    n;
  float  chargefC[10][18];
  int    adc[10][18];
  int    zside[18];
  int    section[18];
  int    channel[18];
};

//
// class declaration
//

class ZDCTreeProducer : public edm::EDAnalyzer {
public:
  explicit ZDCTreeProducer(const edm::ParameterSet&);
  ~ZDCTreeProducer();

  math::XYZPoint getPosition(const DetId &id, reco::Vertex::Point& vtx);
  double getEt(math::XYZPoint& pos, double energy);
  double getEt(const DetId &id, double energy);
  double getEta(const DetId &id);
  double getPhi(const DetId &id);
  double getPerp(const DetId &id);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------

  MyZDCRecHit zdcRecHit;
  MyZDCDigi zdcDigi;

  TNtuple* nt;
  TTree* zdcRecHitTree;
  TTree* zdcDigiTree;

  int nZdcTs_;
  bool calZDCDigi_;

  edm::Service<TFileService> fs;
  edm::ESHandle<CaloGeometry> geo;

  //edm::EDGetTokenT<ZDCDigiCollection> zdcDigiSrc_;
  edm::InputTag zdcDigiSrc_;
  edm::EDGetTokenT<ZDCRecHitCollection> zdcRecHitSrc_;

  bool doZDCRecHit_;
  bool doZDCDigi_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//
constexpr double cone2 = 0.5 * 0.5;

//
// constructors and destructor
//
ZDCTreeProducer::ZDCTreeProducer(const edm::ParameterSet& iConfig) {
  //now do what ever initialization is needed
  doZDCRecHit_ = iConfig.getParameter<bool>("doZDCRecHit");
  doZDCDigi_ = iConfig.getParameter<bool>("doZDCDigi");

  if (doZDCDigi_)
//    zdcDigiSrc_ = consumes<ZDCDigiCollection>(iConfig.getParameter<edm::InputTag>("zdcDigiSrc"));
	 zdcDigiSrc_ = iConfig.getParameter<edm::InputTag>("zdcDigiSrc");
    consumes<QIE10DigiCollection>(zdcDigiSrc_);

  if (doZDCRecHit_)
    zdcRecHitSrc_ = consumes<ZDCRecHitCollection> (iConfig.getParameter<edm::InputTag>("zdcRecHitSrc"));

  nZdcTs_ = iConfig.getParameter<int>("nZdcTs");
  calZDCDigi_ = iConfig.getParameter<bool>("calZDCDigi");
}


ZDCTreeProducer::~ZDCTreeProducer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
ZDCTreeProducer::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
  iSetup.get<CaloGeometryRecord>().get(geo);
  if(doZDCRecHit_){
    edm::Handle<ZDCRecHitCollection> zdcrechits;
    ev.getByToken(zdcRecHitSrc_,zdcrechits);

    int nhits = 0;
    for (auto const& rh : *zdcrechits) {
      HcalZDCDetId zdcid = rh.id();
      if (nhits  < 18) {
        zdcRecHit.e[nhits] = rh.energy();
        zdcRecHit.zside[nhits] = zdcid.zside();
        zdcRecHit.section[nhits] = zdcid.section();
        zdcRecHit.channel[nhits] = zdcid.channel();
        zdcRecHit.saturation[nhits] = static_cast<int>( rh.flagField(HcalCaloFlagLabels::ADCSaturationBit) );
      }

      nhits++;
    } // end loop zdc rechits

    zdcRecHit.n = nhits;
    zdcRecHitTree->Fill();
  }

  if(doZDCDigi_){
//    edm::Handle<ZDCDigiCollection> zdcdigis;
	 edm::Handle<QIE10DigiCollection> zdcdigis;
   ev.getByLabel(zdcDigiSrc_,zdcdigis);

    edm::ESHandle<HcalDbService> conditions;
    iSetup.get<HcalDbRecord>().get(conditions);

    int nhits = 0;
//    for (auto const& rh : *zdcdigis)  {

		std::cout << "zdcdigis->size() : " << zdcdigis->size() << std::endl;
		std::cout << "zdcdigis->samples() : " << zdcdigis->samples() << std::endl;


      for (auto it = zdcdigis->begin(); it != zdcdigis->end(); it++) {
		std::cout << "--- nhits : " << nhits << std::endl;

       const QIE10DataFrame digi = static_cast<const QIE10DataFrame>(*it);


      HcalZDCDetId zdcid = digi.id();


      CaloSamples caldigi;

		 //const ZDCDataFrame & rh = (*zdcdigis)[it];
      if(calZDCDigi_){
	std::cout << "###DIGI1?" << std::endl;


        const HcalQIECoder* qiecoder = conditions->getHcalCoder(zdcid);
        const HcalQIEShape* qieshape = conditions->getHcalShape(qiecoder);
        HcalCoderDb coder(*qiecoder, *qieshape);
//        coder.adc2fC(rh,caldigi);
			 coder.adc2fC(digi,caldigi);
	std::cout << "###--- END DIGI1?" << std::endl;


      }

      if (nhits  < 18) {
        zdcDigi.zside[nhits] = zdcid.zside();
        zdcDigi.section[nhits] = zdcid.section();
        zdcDigi.channel[nhits] = zdcid.channel();

        for (int ts = 0; ts < digi.samples(); ts++) {
	std::cout << "###DIGI2? --- ts : " << ts << std::endl;


          zdcDigi.chargefC[ts][nhits] = calZDCDigi_ ? caldigi[ts] : QWAna::ZDC2018::QIE10_nominal_fC[ digi[ts].adc() ];
          zdcDigi.adc[ts][nhits] = digi[ts].adc();
	std::cout << "###--- END DIGI2?" << std::endl;


        }
      }
      nhits++;
    } // end loop zdc rechits
	std::cout << "###DIGI END?" << std::endl;


    zdcDigi.n = nhits;
    zdcDigiTree->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void
ZDCTreeProducer::beginJob()
{

  if(doZDCRecHit_){
    zdcRecHitTree = fs->make<TTree>("zdcrechit", "zdc");
    zdcRecHitTree->Branch("n",&zdcRecHit.n,"n/I");
    zdcRecHitTree->Branch("e",zdcRecHit.e,"e[n]/F");
    zdcRecHitTree->Branch("saturation",zdcRecHit.saturation,"saturation[n]/F");
    zdcRecHitTree->Branch("zside",zdcRecHit.zside,"zside[n]/I");
    zdcRecHitTree->Branch("section",zdcRecHit.section,"section[n]/I");
    zdcRecHitTree->Branch("channel",zdcRecHit.channel,"channel[n]/I");
  }

  if(doZDCDigi_){
    zdcDigiTree = fs->make<TTree>("zdcdigi", "zdc");
    zdcDigiTree->Branch("n",&zdcDigi.n,"n/I");
    zdcDigiTree->Branch("zside",zdcDigi.zside,"zside[n]/I");
    zdcDigiTree->Branch("section",zdcDigi.section,"section[n]/I");
    zdcDigiTree->Branch("channel",zdcDigi.channel,"channel[n]/I");

    for( int i=0; i<nZdcTs_;i++){
      TString adcTsSt("adcTs"), chargefCTsSt("chargefCTs");
      adcTsSt+=i; chargefCTsSt+=i;

      zdcDigiTree->Branch(adcTsSt,zdcDigi.adc[i],adcTsSt+"[n]/I");
      zdcDigiTree->Branch(chargefCTsSt,zdcDigi.chargefC[i],chargefCTsSt+"[n]/F");
    }
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void
ZDCTreeProducer::endJob() {
}

math::XYZPoint ZDCTreeProducer::getPosition(const DetId &id, reco::Vertex::Point& vtx) {
  const GlobalPoint& pos = geo->getPosition(id);
  return math::XYZPoint(pos.x() - vtx.x(), pos.y() - vtx.y(), pos.z() - vtx.z());
}

inline double ZDCTreeProducer::getEt(math::XYZPoint& pos, double energy) {
  return energy * sin(pos.theta());
}

inline double ZDCTreeProducer::getEt(const DetId &id, double energy) {
  const GlobalPoint& pos = geo->getPosition(id);
  return energy * sin(pos.theta());
}

inline double ZDCTreeProducer::getEta(const DetId &id) {
  const GlobalPoint& pos = geo->getPosition(id);
  return pos.eta();
}

inline double ZDCTreeProducer::getPhi(const DetId &id) {
  const GlobalPoint& pos = geo->getPosition(id);
  return pos.phi();
}

inline double ZDCTreeProducer::getPerp(const DetId &id) {
  const GlobalPoint& pos = geo->getPosition(id);
  return pos.perp();
}


//define this as a plug-in
DEFINE_FWK_MODULE(ZDCTreeProducer);
