/**
 * @file
 * @brief Definition of module EventLoaderVMMSDAT
 *
 * @copyright Copyright (c) 2015-2024 CERN and the Corryvreckan authors.
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
    /** @ingroup Modules
     * @brief Module to do function
     *
     * More detailed explanation of module
     */
    class EventLoaderVMMSDAT : public Module {

    public:
        /**
         * @brief Constructor for this unique module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param detectors Vector of pointers to the detectors
         */
        EventLoaderVMMSDAT(Configuration& config, std::vector<std::shared_ptr<Detector>> detectors);

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
        void loop();
		bool triggered(std::string name, double pos0, double pos1, double charge);
      
    private:
        int eventNumber_;
        std::map<std::string, TH2D*>  hitmap_;
        std::map<int, std::string>    detector_map_;

        std::map<std::string, double> det_pos_scale_;   // should be moved to detector geometry definition class
        std::map<std::string, double> det_pos_shift_x_;   // should be moved to detector geometry definition class
        std::map<std::string, double> det_pos_shift_y_;   // should be moved to detector geometry definition class

        ClusterVector                 runClusters_;

        std::string input_file_name_;
        std::string tree_name_;
        std::string detector_trigger_;
        double channel_trigger_;
        double charge_trigger_;
        
        std::vector<std::string> detector_required_;
        std::unique_ptr<TFile> input_file_;
        TTree       *event_tree_;
        Long64_t number_clusters_to_read_;
        Long64_t number_clusters_to_skip_;
        double time_window_;
        size_t currentCluster_;

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
        unsigned short size0_;
        unsigned short size1_;
        double pos0_;
        double pos1_;
        double pos0_charge2_;
        double pos1_charge2_;
        double pos0_utpc_;
        double pos1_utpc_;
        double pos0_algo_;
        double pos1_algo_;

        std::string position_algorithm_;
        std::string time_algorithm_;
        std::string time_choice_;
        std::string charge_choice_;
        std::string size_choice_;
        bool sort_clusters_;

    // List of branches
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
       TBranch* b_size0;
       TBranch* b_size1;  
       TBranch* b_pos0;
       TBranch* b_pos1;
       TBranch* b_pos0_charge2;
       TBranch* b_pos1_charge2;
       TBranch* b_pos0_utpc;
       TBranch* b_pos1_utpc;
       TBranch* b_pos0_algo;
       TBranch* b_pos1_algo;              
       std::map<std::string, TH2D*> clustermap_{};
       std::map<std::string, TH1D*> numberClustermap_{};
       TH1D*  totalClusters_;
    };

} // namespace corryvreckan
