#pragma once
#include <vector>
struct Hit {
    uint32_t id;
    uint32_t event;
    uint8_t det;
    uint8_t plane;
    uint8_t fec;
    uint8_t vmm;
    double readout_time;
    double time;
    uint8_t ch;
    uint16_t pos;
    uint16_t bcid;
    uint16_t tdc;
    uint16_t adc;
    bool over_threshold;
    float chip_time;
};

struct ClusterPlane {
    uint32_t id;
    uint8_t det;
    uint8_t plane;
    uint16_t size;
    uint16_t adc;
    double time;
    double time_utpc;
    double time_charge2;
    double time_algo;
    double pos;
    double pos_utpc;
    double pos_charge2;
    double pos_algo;
    bool plane_coincidence;
    uint16_t max_delta_time;
    uint16_t max_missing_strip;
    uint16_t span_cluster;
    std::vector<double> strips;
    std::vector<double> times;
    std::vector<double> adcs;

};

struct ClusterDetector {
    uint32_t id;
    uint32_t id0;
    uint32_t id1;
    uint8_t det;
    uint16_t size0;
    uint16_t size1;
    uint16_t adc0;
    uint16_t adc1;
    double pos0;
    double pos1;
    double pos2;
    double time0;
    double time1;
    double pos0_utpc;
    double pos1_utpc;
    double pos2_utpc;
    double time0_utpc;
    double time1_utpc;
    double pos0_charge2;
    double pos1_charge2;
    double pos2_charge2;
    double time0_charge2;
    double time1_charge2;
    double pos0_algo;
    double pos1_algo;
    double pos2_algo;
    double time0_algo;
    double time1_algo;
    double dt0;
    double dt1;
    double delta_plane;
    uint16_t span_cluster0;
    uint16_t span_cluster1;
    uint16_t max_delta_time0;
    uint16_t max_delta_time1;
    uint16_t max_missing_strip0;
    uint16_t max_missing_strip1;
    std::vector<double> strips0;
    std::vector<double> times0;
    std::vector<double> strips1;
    std::vector<double> times1;
    std::vector<double> adcs0;
    std::vector<double> adcs1;
};



struct Track {
    uint32_t id;
    std::vector<uint8_t> det;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> t;

    std::vector<double> x_utpc;
    std::vector<double> y_utpc;
    std::vector<double> z_utpc;
    std::vector<double> t_utpc;

    std::vector<double> x_charge2;
    std::vector<double> y_charge2;
    std::vector<double> z_charge2;
    std::vector<double> t_charge2;
};
using std::string;


using HitTuple = std::tuple<double, uint16_t, uint16_t>;
using ClusterTuple = std::tuple<uint16_t, double, uint16_t>;
using HitContainer = std::vector<HitTuple>;
using ClusterContainer = std::vector<ClusterTuple>;

using ClusterVectorPlane = std::vector<ClusterPlane>;
using ClusterVectorDetector = std::vector<ClusterDetector>;
using TrackVector = std::vector<Track>;
using HitVector = std::vector<Hit>;

