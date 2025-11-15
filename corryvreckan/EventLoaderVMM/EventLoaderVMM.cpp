/**
 * @file
 * @brief Implementation of module EventLoaderVMM
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include "EventLoaderVMM.h"

using namespace corryvreckan;

EventLoaderVMM::EventLoaderVMM(Configuration &config, std::vector<std::shared_ptr<Detector>> detectors) : Module(config, std::move(detectors)) {}

void EventLoaderVMM::initialize() {
  for (auto &detector : get_detectors()) {
    LOG(DEBUG) << "Initialise detector " + detector->getName();
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

  position_algorithm_ = config_.get<std::string>("position_algorithm", "cog");
  time_algorithm_ = config_.get<std::string>("time_algorithm", "cog");
  time_choice_ = config_.get<std::string>("time_choice", "time0");
  charge_choice_ = config_.get<std::string>("charge_choice", "adc0");
  size_choice_ = config_.get<std::string>("size_choice", "size0");

  sort_clusters_ = config_.get<bool>("sort_clusters", true);

  for (auto &detm : detector_map_) {
    LOG(INFO) << "Detector map : [id, name]  " << detm.first << "  " << detm.second;
  }

  LOG(INFO) << "Trigger: " << detector_trigger_;
  LOG(INFO) << "Position algorithm (vmm-sdat): " << position_algorithm_;
  LOG(INFO) << "Time algorithm (vmm-sdat): " << time_algorithm_;
  LOG(INFO) << "Time choice for Corryvreckan cluster: " << time_choice_;
  LOG(INFO) << "Charge choice for Corryvreckan cluster: " << charge_choice_;
  LOG(INFO) << "Sort clusters by time: " << sort_clusters_;

  // Initialise member variables
  eventNumber_ = 0;
  // Start at -1, since in run currentCluster is directly incremented to 0
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

  // Check first existence of these branches for extended functionality
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
  LOG(INFO) << "Set to skip " << number_clusters_to_skip_ << " and then read " << number_clusters_to_read_ << " clusters from the input file.";
  LOG(INFO) << "Use time window of " << time_window_ << " to allocate clusters to event.";

  loop();
  input_file_->Close();
}

StatusCode EventLoaderVMM::run(const std::shared_ptr<Clipboard> &clipboard) {

  std::map<std::string, ClusterVector> deviceClusters;
  ++currentCluster_;
  if (static_cast<size_t>(currentCluster_) < runClusters_.size()) {
    while (!triggered(runClusters_[static_cast<size_t>(currentCluster_)]->getDetectorID(), runClusters_[static_cast<size_t>(currentCluster_)]->row(), runClusters_[static_cast<size_t>(currentCluster_)]->column(), runClusters_[static_cast<size_t>(currentCluster_)]->charge())) {
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

  // cluster_before has to be type int, since it can get negative values
  // if it is size_t, it gets huge positive values leading to a segmentation fault
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
    LOG(TRACE) << "Event " << eventNumber_ << ": has " << detClustItr.second.size() << " clusters for device " << detClustItr.first << "  Current trigger detector cluster: " << currentCluster_;
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

  bool all_clusters = true;
  bool no_clusters = true;
  for (auto &detector : get_detectors()) {
    if ((!detector->hasRole(DetectorRole::AUXILIARY)) && (!detector->hasRole(DetectorRole::DUT))) {
      if (nclsev[detector->getName()] == 0) {
        all_clusters = false;
      }
      if (nclsev[detector->getName()] >= 1) {
        no_clusters = false;
      }
    }
  }

  if (all_clusters) {
    totalClusters_->Fill(-1.0);
  }
  if (no_clusters) {
    totalClusters_->Fill(-2.0);
  }

  // Return value telling analysis to keep running
  return StatusCode::Success;
}

void EventLoaderVMM::finalize(const std::shared_ptr<ReadonlyClipboard> &) {}

void EventLoaderVMM::loop() {
  std::map<std::string, size_t> detposmap;
  size_t ind = 0;
  for (auto &detector : get_detectors()) {
    detposmap[detector->getName()] = ind++;
    Configuration config = detector->getConfiguration();
    if (config.has("number_of_pixels") && config.has("pixel_pitch")) {
      det_size_x_[detector->getName()] = detector->getSize().x();
      det_size_y_[detector->getName()] = detector->getSize().y();
      det_pitch_x_[detector->getName()] = detector->getPitch().x();
      det_pitch_y_[detector->getName()] = detector->getPitch().y();
      LOG(INFO) << "Detector: " << detector->getName() << ": scale_x: " << det_pitch_x_[detector->getName()] << ", shift_x: " <<0.5*det_pitch_x_[detector->getName()] -0.5 * det_size_x_[detector->getName()];
      LOG(INFO) << "Detector: " << detector->getName() << ": scale_y: " << det_pitch_y_[detector->getName()] << ", shift_y: " << 0.5*det_pitch_y_[detector->getName()] -0.5 * det_size_y_[detector->getName()];
    } else {
      if (!detector->hasRole(DetectorRole::AUXILIARY)) {
        LOG(ERROR) << "Detector " << detector->getName() << " does not have number of pixels and pitch defined in geometry!";
      }
      det_size_x_[detector->getName()] = 0.0;
      det_size_y_[detector->getName()] = 0.0;
      det_pitch_x_[detector->getName()] = 0.0;
      det_pitch_y_[detector->getName()] = 0.0;
    }
    if (det_size_x_[detector->getName()] == det_pitch_x_[detector->getName()]) {
      if (det_size_y_[detector->getName()] != det_pitch_y_[detector->getName()]) {
      	LOG(INFO) << "1D-Detector: " << detector->getName() << ": x=0";
      }
      else {
      	LOG(INFO) << "0D-Detector: " << detector->getName() << ": x=0, y=0";
      }
    } else {
      if (det_size_y_[detector->getName()] == det_pitch_y_[detector->getName()]) {
        LOG(INFO) << "1D-Detector: " << detector->getName() << ": y=0";
      }
    }
  }

  Long64_t last_entry = number_clusters_to_skip_ + number_clusters_to_read_;
  if (last_entry > event_tree_->GetEntries())
    last_entry = event_tree_->GetEntries();
  LOG(INFO) << "Reading data from the tree starting from entry " << number_clusters_to_skip_ << " till " << last_entry;

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
    
    if (det_size_x_[detName] == det_pitch_x_[detName]) {
      the_pos0=0;
    }
    if (det_size_y_[detName] == det_pitch_y_[detName]) {
      the_pos1=0;
    }
    
    double scale_x = det_pitch_x_[detName];
    double scale_y = det_pitch_y_[detName];
    //The correct shift also shifts the pixel position to the center of the pixel
    double shift_x = 0.5*det_pitch_x_[detName]-0.5 * det_size_x_[detName];
    double shift_y = 0.5*det_pitch_y_[detName]-0.5 * det_size_y_[detName];

    PositionVector3D<Cartesian3D<double>> positionLocal(the_pos0 * scale_x + shift_x, the_pos1 * scale_y + shift_y, 0.0);
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

bool EventLoaderVMM::triggered(std::string name, double pos0, double pos1, double charge) {
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
