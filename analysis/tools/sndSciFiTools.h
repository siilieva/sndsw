#pragma once

#include "TClonesArray.h"

namespace snd {
  namespace analysis_tools {

    // Function to get the number of hits per station in SciFi
    void getSciFiHitsPerStation(TClonesArray * digiHits, std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits);

    // Function to get the total number of SciFi hits
    int getTotalSciFiHits(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits);
    int getTotalSciFiHits(TClonesArray * digiHits);

    // Function to get the SciFi fractional hits per plane
    std::vector<float> getFractionalHitsPerScifiPlane(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits);
    std::vector<float> getFractionalHitsPerScifiPlane(TClonesArray * digiHits);

    // Function to find the station where the interaction ocurred by checking the station at which the cumulative hit fraction exceeds a threshold
    int findScifiStation(std::vector<int> &horizontal_hits, std::vector<int> &vertical_hits, float threshold);
    int findScifiStation(TClonesArray * digiHits, float threshold);

  }
}
