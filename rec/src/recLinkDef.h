// NA6PCCopyright

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class PID + ;
#pragma link C++ class MatBudget + ;
#pragma link C++ class Propagator + ;
#pragma link C++ class NA6PBaseCluster + ;
#pragma link C++ class std::vector < NA6PBaseCluster> + ;
#pragma link C++ class NA6PTrackPar + ;
#pragma link C++ class NA6PTrackParCov + ;
#pragma link C++ class std::vector < NA6PTrackPar> + ;
#pragma link C++ class std::vector < NA6PTrackParCov> + ;
#pragma link C++ class NA6PVerTelCluster + ;
#pragma link C++ class std::vector < NA6PVerTelCluster> + ;
#pragma link C++ class NA6PMuonSpecCluster + ;
#pragma link C++ class std::vector < NA6PMuonSpecCluster> + ;
#pragma link C++ class NA6PTrack + ;
#pragma link C++ class std::vector < NA6PTrack> + ;
#pragma link C++ class NA6PVertex + ;
#pragma link C++ class std::vector < NA6PVertex> + ;
#pragma link C++ class NA6PMatch + ;
#pragma link C++ class std::vector < NA6PMatch> + ;

#pragma link C++ class NA6PLine + ;
#pragma link C++ class NA6PHelixHelper + ;
#pragma link C++ class Seeder + ;

#pragma link C++ class NA6PReconstruction + ;
#pragma link C++ class NA6PVerTelReconstruction + ;
#pragma link C++ class NA6PMuonSpecReconstruction + ;
#pragma link C++ class NA6PVerTelClusterizer + ;
#pragma link C++ class NA6PMatching + ;
#pragma link C++ class NA6PVertexerTracklets + ;
#pragma link C++ class NA6PVertexerTracks + ;
#pragma link C++ class NA6PTrackerCA + ;
#pragma link C++ class NA6PFastTrackFitter + ;

#pragma link C++ class NA6PRecoParam + ;
#pragma link C++ class na6p::conf::ConfigurableParamHelper < NA6PRecoParam> + ;


#endif
