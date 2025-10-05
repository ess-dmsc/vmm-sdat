/**
 * @file
 * @brief Implementation of module EventLoaderVMMSDAT
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventLoaderVMMSDAT.h"

using namespace corryvreckan;

EventLoaderVMMSDAT::EventLoaderVMMSDAT(Configuration &config, std::vector<std::shared_ptr<Detector>> detectors)
    : Module(config, std::move(detectors)) {}

void EventLoaderVMMSDAT::initialize() {
  for (auto &detector : get_detectors()) {
    //LOG(DEBUG) << "Initialise for detector " + detector->getName();
    std::string detectorID = detector->getName();
    std::string hhname = detectorID + "_ClusterMap";
    std::string title = detectorID + ": clustermap;x [px];y [px];events";
    int npxx = detector->nPixels().X();
    if (!npxx) {
      npxx = 256;
    }
    int npxy = detector->nPixels().Y();
    if (!npxy) {
      npxy = 256;
    }
    clustermap_[detectorID] = new TH2D(hhname.c_str(), title.c_str(), npxx, 0.0, npxx, npxy, 0.0, npxy);

    hhname = detectorID + "_NClustersPerEvent";
    title = detectorID + ": _NClustersPerEvent;N_Clusters;events";
    numberClustermap_[detectorID] = new TH1D(hhname.c_str(), title.c_str(), 20, 0.0, 20.0);
  }
  totalClusters_ = new TH1D("NClusterTotalEvent", "NClusterTotalEvent", 55, -5.0, 50.0);

  detector_map_ = config_.getMap<int, std::string>("detector_map");
  detector_trigger_ = config_.get<std::string>("detector_trigger", "TRIGGER");
  channel_trigger_ = config_.get<double>("channel_trigger", -1.0);
  charge_trigger_ = config_.get<double>("charge_trigger", -1.0);
  bool triggerDetectorFound = false;
  for (auto &det : detector_map_) {
    if (det.second == detector_trigger_) {
      triggerDetectorFound = true;
      break;
    }
  }
  if (triggerDetectorFound == false) {
    LOG(ERROR) << detector_trigger_ << " is not in the list of detectors!";
    return;
  }

  if (!config_.has("detector_required")) {
    LOG(ERROR) << "Required detectors have to be defined!";
    return;
  }
  detector_required_ = config_.getArray<std::string>("detector_required");

  position_algorithm_ = config_.get<std::string>("position_algorithm", "cog");
  time_algorithm_ = config_.get<std::string>("time_algorithm", "cog");
  time_choice_ = config_.get<std::string>("time_choice", "time0");
  charge_choice_ = config_.get<std::string>("charge_choice", "adc0");
  size_choice_ = config_.get<std::string>("size_choice", "size0");

  sort_clusters_ = config_.get<bool>("sort_clusters", true);

  if (config_.has("detector_position_scale_map")) {
    det_pos_scale_ = config_.getMap<std::string, double>("detector_position_scale_map");
  }
  if (config_.has("detector_position_shift_x_map")) {
    det_pos_shift_x_ = config_.getMap<std::string, double>("detector_position_shift_x_map");
  }
  if (config_.has("detector_position_shift_y_map")) {
    det_pos_shift_y_ = config_.getMap<std::string, double>("detector_position_shift_y_map");
  }

  for (auto &detm : detector_map_) {
    if (det_pos_scale_.find(detm.second) == det_pos_scale_.end()) {
      LOG(WARNING) << "Position scale factor for detector " << detm.second << " is not set in config file, using 1.0";
      det_pos_scale_[detm.second] = 1.0;
    }

    if (det_pos_shift_x_.find(detm.second) == det_pos_shift_x_.end()) {
      LOG(WARNING) << "Position shift in x for detector " << detm.second << " is not set in config file, using 0.0";
      det_pos_shift_x_[detm.second] = 0.0;
    }

    if (det_pos_shift_y_.find(detm.second) == det_pos_shift_y_.end()) {
      LOG(WARNING) << "Position shift in y for detector " << detm.second << " is not set in config file, using 0.0";
      det_pos_shift_y_[detm.second] = 0.0;
    }
  }

  for (auto &detm : detector_map_) {
    LOG(INFO) << "Detector map : [id, name]  " << detm.first << "  " << detm.second;
  }
  for (auto &detm : det_pos_scale_) {
    LOG(INFO) << "Detector position scale factor  : [name, scale]  " << detm.first << "  " << detm.second;
  }
  for (auto &detm : det_pos_shift_x_) {
    LOG(INFO) << "Detector position shift x : [name, shift_x]  " << detm.first << "  " << detm.second;
  }
  for (auto &detm : det_pos_shift_y_) {
    LOG(INFO) << "Detector position shift y : [name, shift_y]  " << detm.first << "  " << detm.second;
  }
  for (auto &detv : detector_required_) {
    LOG(INFO) << "Required detector  " << detv;
  }
  LOG(INFO) << "Trigger: " << detector_trigger_;
  LOG(INFO) << "Position algorithm (vmm-sdat): " << position_algorithm_;
  LOG(INFO) << "Time algorithm (vmm-sdat): " << time_algorithm_;
  LOG(INFO) << "Time choice for Corryvreckan cluster: " << time_choice_;
  LOG(INFO) << "Charge choice for Corryvreckan cluster: " << charge_choice_;
  LOG(INFO) << "Sort clusters by time: " << sort_clusters_;

  // Initialise member variables
  eventNumber_ = 0;
  //Start at -1, since in run currentCluster is directly incremented to 0
  currentCluster_ = -1;

  input_file_name_ = config_.get<std::string>("input_file");
  tree_name_ = config_.get<std::string>("tree_name", "clusters_detector");
  number_clusters_to_read_ = config_.get<Long64_t>("number_clusters_to_read", -1);
  number_clusters_to_skip_ = config_.get<Long64_t>("number_clusters_to_skip", 0);

  time_window_ = config_.get<double>("time_window", 2000.0);

  input_file_ = std::make_unique<TFile>(input_file_name_.c_str(), "READ");
  if (!input_file_) {
    LOG(ERROR) << "Cannot open input file " << input_file_name_;
    return;
  }

  event_tree_ = input_file_->Get<TTree>(tree_name_.c_str());

  if (!event_tree_) {
    LOG(ERROR) << "Cannot get tree " << tree_name_ << " fdom the inpu file " << input_file_name_;
    return;
  }

  event_tree_->SetBranchAddress("time0", &time0_, &b_time0);
  event_tree_->SetBranchAddress("det", &det_, &b_det);
  event_tree_->SetBranchAddress("adc0", &adc0_, &b_adc0);
  event_tree_->SetBranchAddress("size0", &size0_, &b_size0);
  event_tree_->SetBranchAddress("pos0", &pos0_, &b_pos0);
  event_tree_->SetBranchAddress("pos1", &pos1_, &b_pos1);
  
  //Check first existence of these branches for extended functionality
  if (event_tree_->GetBranch("size1")) {
    event_tree_->SetBranchAddress("size1", &size1_, &b_size1);
  } 
  if (event_tree_->GetBranch("adc1")) {
    event_tree_->SetBranchAddress("adc1", &adc1_, &b_adc1);
  } 
  if (event_tree_->GetBranch("time1")) {
    event_tree_->SetBranchAddress("time1", &time1_, &b_time1);
  }
  if (event_tree_->GetBranch("time0_charge2")) {
    event_tree_->SetBranchAddress("time0_charge2", &time0_charge2_, &b_time0_charge2);
  } 
  if (event_tree_->GetBranch("time1_charge2")) {
    event_tree_->SetBranchAddress("time1_charge2", &time1_charge2_, &b_time1_charge2);
  } 
  if (event_tree_->GetBranch("time0_utpc")) {
    event_tree_->SetBranchAddress("time0_utpc", &time0_utpc_, &b_time0_utpc);
  } 
  if (event_tree_->GetBranch("time1_utpc")) {
    event_tree_->SetBranchAddress("time1_utpc", &time1_utpc_, &b_time1_utpc);
  }
  if (event_tree_->GetBranch("time0_algo")) {
    event_tree_->SetBranchAddress("time0_algo", &time0_algo_, &b_time0_algo);
  } 
  if (event_tree_->GetBranch("time1_algo")) {
    event_tree_->SetBranchAddress("time1_algo", &time1_algo_, &b_time1_algo);
  } 
  if (event_tree_->GetBranch("pos0_charge2")) {
    event_tree_->SetBranchAddress("pos0_charge2", &pos0_charge2_, &b_pos0_charge2);
  }
  if (event_tree_->GetBranch("pos1_charge2")) {
    event_tree_->SetBranchAddress("pos1_charge2", &pos1_charge2_, &b_pos1_charge2);
  }
  if (event_tree_->GetBranch("pos0_utpc")) {
    event_tree_->SetBranchAddress("pos0_utpc", &pos0_utpc_, &b_pos0_utpc);
  }
  if (event_tree_->GetBranch("pos1_utpc")) {
    event_tree_->SetBranchAddress("pos1_utpc", &pos1_utpc_, &b_pos1_utpc);
  }
  if (event_tree_->GetBranch("pos0_algo")) {
    event_tree_->SetBranchAddress("pos0_algo", &pos0_algo_, &b_pos0_algo);
  }
  if (event_tree_->GetBranch("pos1_algo")) {
    event_tree_->SetBranchAddress("pos1_algo", &pos1_algo_, &b_pos1_algo);
  }  
     

  if ((number_clusters_to_read_ < 0) || (number_clusters_to_read_ > event_tree_->GetEntries())) {
    number_clusters_to_read_ = event_tree_->GetEntries();
  }
  LOG(INFO) << "Set to skip " << number_clusters_to_skip_ << " and then read " << number_clusters_to_read_ << " clusters form the input file.";
  LOG(INFO) << "Use time window of " << time_window_ << " for collecting clustres to event.";

  loop();
  input_file_->Close();
}

StatusCode EventLoaderVMMSDAT::run(const std::shared_ptr<Clipboard> &clipboard) {


  std::map<std::string, ClusterVector> deviceClusters;
  ++currentCluster_;
  if (static_cast<size_t>(currentCluster_) < runClusters_.size()) {
    while (!triggered(runClusters_[static_cast<size_t>(currentCluster_)]->getDetectorID(), runClusters_[static_cast<size_t>(currentCluster_)]->row(), runClusters_[static_cast<size_t>(currentCluster_)]->column(),
                      runClusters_[static_cast<size_t>(currentCluster_)]->charge())) {                
      if (static_cast<size_t>(++currentCluster_) == runClusters_.size()) {
        LOG(INFO) << "All data read from the file! EndRun.";
        return StatusCode::Failure;
      }
    }
  } else {
    LOG(INFO) << "All data read from the file! EndRun.";
    return StatusCode::Failure;
  }  
  
  std::shared_ptr<Cluster> ccluster = runClusters_[static_cast<size_t>(currentCluster_)];
  double time_window = time_window_;
    
  // These two loops will add all clusters in the time window around selected one, assuming the vector is time sorted
  std::map<std::string, size_t> nclsev;
  std::for_each(detector_map_.begin(), detector_map_.end(), [&](const auto &pp) { nclsev[pp.second] = 0; });

  double eve_timestamp = ccluster->timestamp();
  size_t cluster_after = static_cast<size_t>(currentCluster_);
  while (cluster_after < runClusters_.size() && std::abs(runClusters_[cluster_after]->timestamp() - eve_timestamp) <= time_window) {
    std::shared_ptr<Cluster> iclust = runClusters_[cluster_after];
    deviceClusters[iclust->getDetectorID()].push_back(iclust);
    ++nclsev[iclust->getDetectorID()];
    ++cluster_after;
  }

  //cluster_before has to be type int, since it can get negative values
  //if it is size_t, it gets huge positive values leading to a segmentation fault
  int cluster_before = currentCluster_ - 1;
  while (cluster_before >= 0 && (std::abs(runClusters_[static_cast<size_t>(cluster_before)]->timestamp() - eve_timestamp) <= time_window)) {
    std::shared_ptr<Cluster> jclust = runClusters_[static_cast<size_t>(static_cast<size_t>(cluster_before))];
    deviceClusters[jclust->getDetectorID()].push_back(jclust);
    ++nclsev[jclust->getDetectorID()];
    --cluster_before;
  }

  for (const auto &detClustItr : deviceClusters) {
    // Add dummy pixels as they are needed by the other modules of Corryvreckan
    PixelVector pixels;
    for (auto &dcluster : detClustItr.second) {
      int cl_size = static_cast<int>(std::round(dcluster->widthX()));
      int cl_col = static_cast<int>(std::round(dcluster->column()));
      int cl_row = static_cast<int>(std::round(dcluster->row()));
      double cl_charge = dcluster->charge();
      cl_size = 1; // probably one dummy pixel is enough for algorithms to work.
      for (int pi = 0; pi < cl_size; ++pi) {
        auto pixel = std::make_shared<Pixel>(std::string(detClustItr.first), cl_col, cl_row, cl_charge, cl_charge, eve_timestamp);
        pixels.push_back(pixel);
        dcluster->addPixel(pixel.get());
      }
    }
    clipboard->putData(pixels, detClustItr.first);
    clipboard->putData(detClustItr.second, detClustItr.first);
    LOG(INFO) << "Event " << eventNumber_ << ": has " << detClustItr.second.size() << " clusters for device " << detClustItr.first
              << "  Current trigger detector cluster: " << currentCluster_;
  }
  // Increment event counter
  eventNumber_++;

  size_t nclstot = 0;
  clipboard->putEvent(std::make_shared<Event>(eve_timestamp - time_window, eve_timestamp + time_window));
  // Fill the histograms with the number of clusters in each detector, total and other selections
  std::for_each(nclsev.begin(), nclsev.end(), [&](const auto &pp) {
    numberClustermap_[pp.first]->Fill(static_cast<double>(pp.second));
    nclstot += pp.second;
  });
  totalClusters_->Fill(static_cast<double>(nclstot));

  bool all_zero = true;
  bool all_at_least_one = true;
  for (size_t n = 0; n < detector_required_.size(); n++) {
    if (nclsev[detector_required_[n]] != 0) {
      all_zero = false;
    }
    if (nclsev[detector_required_[n]] < 1) {
      all_at_least_one = false;
    }
  }
  if (all_zero) {
    totalClusters_->Fill(-1.0);
  }
  if (all_at_least_one) {
    totalClusters_->Fill(-2.0);
  }

  // Return value telling analysis to keep running
  return StatusCode::Success;
}

void EventLoaderVMMSDAT::finalize(const std::shared_ptr<ReadonlyClipboard> &) { 
//LOG(DEBUG) << "Analysed " << eventNumber_ << " events"; 
}

void EventLoaderVMMSDAT::loop() {
  std::map<std::string, size_t> detposmap;
  size_t ind = 0;
  for (auto &detector : get_detectors()) {
    detposmap[detector->getName()] = ind++;
    Configuration config = detector->getConfiguration();
    double size_x = 0;
    double size_y = 0;
    double pitch_x = 0;
    if (config.has("number_of_pixels") && config.has("pixel_pitch")) {
      size_x = detector->getSize().x();
      size_y = detector->getSize().y();
      det_pos_shift_x_[detector->getName()] = (-0.5 * size_x);
      det_pos_shift_y_[detector->getName()] = (-0.5 * size_y);
      pitch_x = detector->getPitch().x();
      det_pos_scale_[detector->getName()] = pitch_x;
    } else {
      det_pos_scale_[detector->getName()] = 1.0;
      det_pos_shift_x_[detector->getName()] = 0.0;
      det_pos_shift_y_[detector->getName()] = 0.0;
    }
  }

  Long64_t last_entry = number_clusters_to_skip_ + number_clusters_to_read_;
  if (last_entry > event_tree_->GetEntries())
    last_entry = event_tree_->GetEntries();
  LOG(INFO) << "Reading data from the tree starting form entry " << number_clusters_to_skip_ << " till " << last_entry;

  for (Long64_t i = number_clusters_to_skip_; i < last_entry; i++) {
    event_tree_->GetEntry(i);
    if (detector_map_.find(det_) == detector_map_.end())
      continue;

    auto cluster = std::make_shared<Cluster>();
    std::string detName = detector_map_[det_];

    double the_charge = 0;
    double the_pos0 = 0;
    double the_pos1 = 0;
    double the_time = 0;
    double the_time0 = 0;
    double the_time1 = 0;

    if (position_algorithm_ == "charge2") {
      the_pos0 = pos0_charge2_;
      the_pos1 = pos1_charge2_;
    } else if (position_algorithm_ == "algo") {
      the_pos0 = pos0_algo_;
      the_pos1 = pos1_algo_;
    } else if (position_algorithm_ == "utpc") {
      the_pos0 = pos0_utpc_;
      the_pos1 = pos1_utpc_;
    } else {
      the_pos0 = pos0_;
      the_pos1 = pos1_;
    }

    if (time_algorithm_ == "charge2") {
      the_time0 = time0_charge2_;
      the_time1 = time1_charge2_;
    } else if (time_algorithm_ == "algo") {
      the_time0 = time0_algo_;
      the_time1 = time1_algo_;
    } else if (time_algorithm_ == "utpc") {
      the_time0 = time0_utpc_;
      the_time1 = time1_utpc_;
    } else {
      the_time0 = time0_;
      the_time1 = time1_;
    }

    if (time_choice_ == "plane1") {
      the_time = the_time1;
    } else if (time_choice_ == "both") {
      the_time = 0.5 * (the_time0 + the_time1);
    } else {
      the_time = the_time0;
    }

    if (charge_choice_ == "plane1") {
      the_charge = static_cast<double>(adc1_);
    } else if (charge_choice_ == "both") {
      the_charge = 0.5 * (adc0_ + adc1_);
    } else {
      the_charge = static_cast<double>(adc0_);
    }

    if (size_choice_ == "plane1") {
      cluster->setWidth(1.0, size1_);
    } else if (size_choice_ == "both") {
      cluster->setWidth(size0_, size1_);
    } else {
      cluster->setWidth(size0_, 1.0);
    }

    cluster->setTimestamp(the_time);
    cluster->setDetectorID(detName);
    cluster->setCharge(the_charge);

    PositionVector3D<Cartesian3D<double>> positionLocal(pos0_ * det_pos_scale_[detName] + det_pos_shift_x_[detName],
                                                        pos1_ * det_pos_scale_[detName] + det_pos_shift_y_[detName], 0.0);
    std::shared_ptr<Detector> m_detector = get_detectors().at(detposmap[detName]);
    auto positionGlobal = m_detector->localToGlobal(positionLocal);
    LOG(TRACE) << "Detector: " << detName << ": cluster local position: " << positionLocal;
    LOG(TRACE) << "Detector: " << detName << ": cluster global position: " << positionGlobal;
    clustermap_[detName]->Fill(the_pos0, the_pos1);
    cluster->setClusterCentreLocal(positionLocal);
    cluster->setClusterCentre(positionGlobal);

    cluster->setColumn(the_pos0);
    cluster->setRow(the_pos1);
    cluster->setError(m_detector->getSpatialResolution(the_pos0, the_pos1));

    TMatrixD errorMatrix(3, 3);
    errorMatrix(0, 0) = 1.0;
    errorMatrix(1, 1) = 1.0;
    errorMatrix(2, 2) = 1.0;
    cluster->setErrorMatrixGlobal(errorMatrix);

    runClusters_.push_back(cluster);
  }

  // Sort clusters by time
  if (sort_clusters_) {
    LOG(INFO) << "Data is read. Start sorting.";
    auto compareByTime = [&](const std::shared_ptr<Cluster> a, const std::shared_ptr<Cluster> b) -> bool { return a->timestamp() < b->timestamp(); };
    std::sort(runClusters_.begin(), runClusters_.end(), compareByTime);
    LOG(INFO) << "Data should be sorted.";
  } else {
    LOG(INFO) << "Data is read.";
  }
}

bool EventLoaderVMMSDAT::triggered(std::string name, double pos0, double pos1, double charge) {
  if (name != detector_trigger_) {
    return false;
  }
  if (channel_trigger_ != -1.0) {
    if (channel_trigger_ < pos0 - 1) {
      return false;
    }
    if (channel_trigger_ > pos0 + 1) {
      return false;
    }
    if (channel_trigger_ < pos1 - 1) {
      return false;
    }
    if (channel_trigger_ > pos1 + 1) {
      return false;
    }
  }
  if (charge_trigger_ != -1.0) {
    if (charge_trigger_ < charge) {
      return false;
    }
  }
  return true;
}
