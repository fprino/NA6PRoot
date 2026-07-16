#include <cstdint>
#include "NA6PEventReader.h"

void readEvent(
    const char* fileNameVerTel         = "VerticesVerTel.root",
    const char* fileNameTracksVerTel   = "TracksVerTel.root",
    const char* fileNameTracksMuonSpec = "TracksMuonSpec.root",
    const char* fileNameTracksMatching = "TracksMatching.root",
    const char* fileNameMC             = "MCKine.root",
    bool readVerTel         = true,
    bool readTracksVerTel   = true,
    bool readTracksMuonSpec = true,
    bool readTracksMatching = true,
    bool readMC             = true
)
{
  NA6PEventReader reader(fileNameVerTel, fileNameTracksVerTel, fileNameTracksMuonSpec, fileNameTracksMatching, fileNameMC, readVerTel, readTracksVerTel, readTracksMuonSpec, readTracksMatching, readMC);

  for (std::int64_t iev = 0; iev < reader.entries(); ++iev) {
    reader.loadEvent(iev);

    const auto& vertices = reader.verticesVerTel();
    const auto& vtTracks = reader.tracksVerTel();
    const auto& msTracks = reader.tracksMuonSpec();
    const auto& matches  = reader.matches();
    const auto& mcParts  = reader.mcParticles();
    const auto& mcHeader = reader.mcHeader();

    std::cout << "event " << iev << " nVtx=" << vertices.size() << " nVT=" << vtTracks.size() << " nMS=" << msTracks.size() << " nMatch=" << matches.size() << " nMC=" << mcParts.size() << " mcVZ=" << mcHeader->getVZ() << std::endl;
  }
}
