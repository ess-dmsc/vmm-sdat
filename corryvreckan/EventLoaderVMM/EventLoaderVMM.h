/**
 * @file
 * @brief Definition of module EventLoaderVMM
 *
 * @copyright Copyright (c) 2015-2025 CERN and the Corryvreckan authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 * SPDX-License-Identifier: MIT
 */

#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <iostream>
#include "core/module/Module.hpp"
#include "objects/Cluster.hpp"
#include "objects/Pixel.hpp"
#include "objects/Track.hpp"

namespace corryvreckan {
  	/**
     * @brief EventLoaderVMM loads clustered VMM3a data from the root tree
     * clusters_detector, that is the output of the vmm-sdat software.
     * https://github.com/ess-dmsc/vmm-sdat
     * 
     * EventLoaderVMM carries out the following steps:
     * 1) Sorts the tree clusters_detector by cluster time
     * 2) Loops over all vmm-sdat clusters
     * 3) Creates corryvreckan clusters using the x/y position and cluster time 
     
     * 3) first sorts all clusters in the root tree by time.
     * and then creates corryvreckan 
     */
    class EventLoaderVMM : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detectors Vector of pointers to the detectors
         */
        EventLoaderVMM(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors);

        /**
         * @brief [Initialise this module]
         */
        void initialize() override;

        /**
         * @brief [Run the function of this module]
         */
        StatusCode run(const std::shared_ptr<Clipboard>& clipboard) override;

        /**
         * @brief [Finalise module]
         */
        void finalize(const std::shared_ptr<ReadonlyClipboard>& clipboard) override;

    private:
    	/**
         * @brief [Loop over all clusters in tree]
         * The function loops over all clusters in the vmm-sdat root tree
         * It reads clusters from the vmm-sdat clusters, create a std::make_shared<Cluster>
         * and puts them into the ClusterVector runClusters_.
         * It then optionally sorts the clusters in runClusters_ by time.
         */
        void loop();
        
        /**
         * @brief [triggered]
         */
		bool triggered(std::string name, double pos0, double pos1, double charge);
      
    private:
        int eventNumber_;
        std::map<std::string, TH2D*>  hitmap_;
        std::map<int, std::string>    detector_map_;
	  
		//Geometry parameters are read from geometry file
        std::map<std::string, double> det_size_x_;
        std::map<std::string, double> det_size_y_;
        std::map<std::string, double> det_pitch_x_;
        std::map<std::string, double> det_pitch_y_;
        std::map<std::string, std::string>    detector_dimension_;
  
        ClusterVector                 runClusters_;

        std::string input_file_name_;
        std::string tree_name_detector_;
        std::string tree_name_plane_;
        std::string detector_trigger_;
        double channel_trigger_;
        double charge_trigger_;
        
        std::unique_ptr<TFile> input_file_;
        TTree       *event_tree_detector_;
        TTree       *event_tree_plane_;
        Long64_t number_clusters_to_read_;
        Long64_t number_clusters_to_skip_;
        double time_window_;
        int currentCluster_;

        double time0_;
        double time1_;
        double time0_charge2_;
        double time1_charge2_;
        double time0_utpc_;
        double time1_utpc_;
        double time0_algo_;
        double time1_algo_;

        unsigned char det_;
        unsigned short adc0_;
        unsigned short adc1_;
        double pos0_;
        double pos1_;
        double pos0_charge2_;
        double pos1_charge2_;
        double pos0_utpc_;
        double pos1_utpc_;
        double pos0_algo_;
        double pos1_algo_;
        
        
        double time_;
        double time_charge2_;
        double time_utpc_;
        double time_algo_;

        unsigned char plane_;
        unsigned short adc_;
        double pos_;
        double pos_charge2_;
        double pos_utpc_;
        double pos_algo_;
        
        std::string position_algorithm_;
        std::string time_algorithm_;
        std::string time_choice_;
        std::string charge_choice_;
        bool sort_clusters_;

        // List of branches tree detector
       TBranch* b_time0;
       TBranch* b_time1;
       TBranch* b_time0_charge2;
       TBranch* b_time1_charge2;
       TBranch* b_time0_utpc;
       TBranch* b_time1_utpc;
       TBranch* b_time0_algo;
       TBranch* b_time1_algo;
       TBranch* b_det;
       TBranch* b_adc0;
       TBranch* b_adc1;    
       TBranch* b_pos0;
       TBranch* b_pos1;
       TBranch* b_pos0_charge2;
       TBranch* b_pos1_charge2;
       TBranch* b_pos0_utpc;
       TBranch* b_pos1_utpc;
       TBranch* b_pos0_algo;
       TBranch* b_pos1_algo;    
       
       // List of branches tree plane
   		TBranch* b_det_p;
  		TBranch* b_plane_p;
  		TBranch* b_adc_p;
  		TBranch* b_time_p;
  		TBranch* b_time_utpc_p;
  		TBranch* b_time_charge2_p;
  		TBranch* b_time_algo_p;
  		TBranch* b_pos_p;
  		TBranch* b_pos_utpc_p;
  		TBranch* b_pos_charge2_p;
  		TBranch* b_pos_algo_p;
                 
       std::map<std::string, TH2D*> clustermap_{};
       std::map<std::string, TH1D*> numberClustermap_{};
       TH1D*  totalClusters_;
    };

} // namespace corryvreckan
