/* Copyright (C) 2018 European Spallation Source, ERIC. See LICENSE file */
//===----------------------------------------------------------------------===//
///
/// \file
///
/// \brief Class for handling calibration files for VMM asics
///
//===----------------------------------------------------------------------===//

//#include <parser/Log.h>
#include <fstream>
#include <parser/CalibrationFile.h>
#include <parser/Trace.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <nlohmann/json.hpp>
#pragma GCC diagnostic pop

using json = nlohmann::json;

//#undef TRC_LEVEL
//#define TRC_LEVEL TRC_L_DEB

namespace Gem {
/// \brief load calibration from file
CalibrationFile::CalibrationFile(std::string jsonfile) : CalibrationFile() {

  if (jsonfile.empty()) {
    return;
  }

  XTRACE(INIT, DEB, "Loading calibration file %s\n", jsonfile.c_str());

  std::ifstream t(jsonfile);
  std::string Jsonstring((std::istreambuf_iterator<char>(t)),
                         std::istreambuf_iterator<char>());
  if (!t.good()) {
    XTRACE(INIT, ERR, "Invalid Json file %s\n", jsonfile.c_str());
    throw std::runtime_error(
        "CalibrationFile error - requested file unavailable.");
  }

  loadCalibration(Jsonstring);
}

/// \brief parse json string with calibration data
void CalibrationFile::loadCalibration(std::string jsonstring) {
  nlohmann::json Root;
  try {
    Root = nlohmann::json::parse(jsonstring);
  } catch (...) {
    XTRACE(INIT, DEB, "Invalid Json: %s\n", jsonstring.c_str());
    throw std::runtime_error("Invalid Json in calibration file.");
  }

  try {
    auto VmmCals = Root["vmm_calibration"];
    for (auto &vmmcal : VmmCals) {
      auto fecid = vmmcal["fecID"].get<size_t>();
      auto vmmid = vmmcal["vmmID"].get<size_t>();
      auto adc_offsets = vmmcal["adc_offsets"];
      auto adc_slopes = vmmcal["adc_slopes"];
      auto time_offsets = vmmcal["time_offsets"];
      auto time_slopes = vmmcal["time_slopes"];
      auto timewalk_as = vmmcal["timewalk_a"];
      auto timewalk_bs = vmmcal["timewalk_b"];
      auto timewalk_cs = vmmcal["timewalk_c"];
      auto timewalk_ds = vmmcal["timewalk_d"];

      XTRACE(INIT, DEB,
             "fecid: %lu, vmmid: %lu, adc_offsets(%lu), adc_slopes(%lu), "
             "time_offsets(%lu), time_slopes(%lu),"
             "timewalk_as(%lu),timewalk_bs(%lu),"
             "timewalk_cs(%lu),timewalk_ds(%lu)\n",
             fecid, vmmid, adc_offsets.size(), adc_slopes.size(),
             time_offsets.size(), time_slopes.size(), timewalk_as.size(),
             timewalk_bs.size(), timewalk_cs.size(), timewalk_ds.size());

      if ((adc_offsets.size() != MAX_CH) or (adc_slopes.size() != MAX_CH) or
          (time_offsets.size() != MAX_CH) or (time_slopes.size() != MAX_CH) or
          (timewalk_as.size() != MAX_CH) or (timewalk_bs.size() != MAX_CH) or
          (timewalk_cs.size() != MAX_CH) or (timewalk_ds.size() != MAX_CH)) {
        XTRACE(INIT, DEB,
               "Invalid channel configuration, skipping for fec {%lu} and vmm "
               "{%lu}",
               fecid, vmmid);
        throw std::runtime_error("Invalid array lengths in calibration file.");
      }

      for (size_t j = 0; j < MAX_CH; j++) {
        addCalibration(fecid, vmmid, j, adc_offsets[j].get<float>(),
                       adc_slopes[j].get<float>(), time_offsets[j].get<float>(),
                       time_slopes[j].get<float>(), timewalk_as[j].get<float>(),
                       timewalk_bs[j].get<float>(), timewalk_cs[j].get<float>(),
                       timewalk_ds[j].get<float>());
      }
    }
  } catch (const std::exception &exc) {
    XTRACE(INIT, DEB, "JSON config - invalid json: {%s}", exc.what());
    throw std::runtime_error("Invalid json in calibration file field.");
  }
}

void CalibrationFile::addCalibration(size_t fecId, size_t vmmId, size_t chNo,
                                     float adc_offset, float adc_slope,
                                     float time_offset, float time_slope,
                                     float timewalk_a, float timewalk_b,
                                     float timewalk_c, float timewalk_d) {
  if (fecId >= Calibrations.size())
    Calibrations.resize(fecId + 1);
  auto &fec = Calibrations[fecId];
  if (vmmId >= fec.size())
    fec.resize(vmmId + 1);
  auto &vmm = fec[vmmId];
  if (chNo >= vmm.size())
    vmm.resize(chNo + 1);

  vmm[chNo] = {adc_offset, adc_slope,  time_offset, time_slope,
               timewalk_a, timewalk_b, timewalk_c,  timewalk_d};
}

const Calibration &CalibrationFile::getCalibration(size_t fecId, size_t vmmId,
                                                   size_t chNo) const {
  if (fecId >= Calibrations.size())
    return NoCorr;
  const auto &fec = Calibrations[fecId];
  if (vmmId >= fec.size())
    return NoCorr;
  const auto &vmm = fec[vmmId];
  if (chNo >= vmm.size())
    return NoCorr;
  return vmm[chNo];
}

std::string CalibrationFile::debug() const {
  std::string ret;
  for (size_t fecID = 0; fecID < Calibrations.size(); ++fecID) {
    const auto &fec = Calibrations[fecID];
    if (!fec.empty()) {
      ret += fmt::format("\n  FEC={}", fecID);
    }
    for (size_t vmmID = 0; vmmID < fec.size(); ++vmmID) {
      const auto &vmm = fec[vmmID];
      if (!vmm.empty()) {
        ret += fmt::format("\n{:>8}{:<10}", "vmm=", vmmID);
      }
      for (size_t chipNo = 0; chipNo < vmm.size(); ++chipNo) {
        const auto &cal = vmm[chipNo];
        if ((chipNo % 8) == 0)
          ret += fmt::format("{:<7}", "\n");
        ret += fmt::format("{:<5}", fmt::format("[{}]", chipNo)) +
               fmt::format("{:>7}", cal.adc_offset) +
               ((cal.adc_slope >= 0.0) ? " +" : " -") +
               fmt::format("{:>5}x    ", std::fabs(cal.adc_slope));
      }
    }
  }
  return ret;
}
} // namespace Gem
