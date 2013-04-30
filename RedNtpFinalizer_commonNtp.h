#ifndef RedNtpFinalizer_commonNtp_h
#define RedNtpFinalizer_commonNtp_h

#include <vector>
#include "TChain.h"
#include "TH1F.h"
#include "TFile.h"
#include "AnalysisJet.h"
#include "RedNtpFinalizer.h"

class RedNtpFinalizer_commonNtp {

 public:

  RedNtpFinalizer_commonNtp( const std::string& analyzerType, const std::string& dataset, const std::string& flags="" );
  virtual ~RedNtpFinalizer_commonNtp();

  void createOutputFile( const std::string& additionalFlags="" );
  virtual void addFile(const std::string& dataseti, const std::string& selection="");

  TChain* get_tree() { return tree_; };
  TFile* get_outFile() { return outFile_; };
  bool get_DEBUG() { return DEBUG_; };

  void clear();

  void set_outFile( const std::string& fileName="", const std::string& suffix="" );
  void set_outputDir( const std::string& outputDir ) { outputDir_ =  outputDir; };
  void set_dataset( const std::string& dataset ) { dataset_ = dataset; };
  void set_redNtpDir( const std::string& redNtpDir ) { redNtpDir_ = redNtpDir; };
  void set_flags( const std::string& flags ) { flags_ = flags; };
  void set_DEBUG( bool DEBUG ) { DEBUG_ = DEBUG; };

  virtual void finalize() = 0;

  TChain* tree_;

  std::string analyzerType_;
  std::string outputDir_;
  std::string redNtpDir_;
  std::string dataset_;
  std::string flags_;

  TFile* outFile_;

  bool DEBUG_;

 private:

};


#endif
