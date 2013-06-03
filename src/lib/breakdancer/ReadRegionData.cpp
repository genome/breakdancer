#include "ReadRegionData.hpp"

ReadRegionData::~ReadRegionData() {
    for (size_t i = 0; i < _regions.size(); ++i) {
        delete _regions[i];
    }
}
