#include "FDSelectionUtils.h"


/*
int FDSelectionUtils::TrueParticleID(const art::Ptr<recob::Hit> hit, bool rollup_unsaved_ids) {
  std::map<int,double> id_to_energy_map;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::vector<sim::TrackIDE> track_ides = bt_serv->HitToTrackIDEs(hit);
  for (unsigned int idIt = 0; idIt < track_ides.size(); ++idIt) {
    int id = track_ides.at(idIt).trackID;
    if (rollup_unsaved_ids) id = std::abs(id);
    double energy = track_ides.at(idIt).energy;
    id_to_energy_map[id]+=energy;
  }
  //Now loop over the map to find the maximum contributor
  double likely_particle_contrib_energy = -99999;
  int likely_track_id = 0;
  for (std::map<int,double>::iterator mapIt = id_to_energy_map.begin(); mapIt != id_to_energy_map.end(); mapIt++){
    double particle_contrib_energy = mapIt->second;
    if (particle_contrib_energy > likely_particle_contrib_energy){
      likely_particle_contrib_energy = particle_contrib_energy;
      likely_track_id = mapIt->first;
    }
  }
  return likely_track_id;
}



int FDSelectionUtils::TrueParticleIDFromTotalTrueEnergy(const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::map<int,double> trackIDToEDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      int id = trackIDs[idIt].trackID;
      if (rollup_unsaved_ids) id = std::abs(id);
      trackIDToEDepMap[id] += trackIDs[idIt].energy;
    }
  }

  //Loop over the map and find the track which contributes the highest energy to the hit vector
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator mapIt = trackIDToEDepMap.begin(); mapIt != trackIDToEDepMap.end(); mapIt++){
    double energy = mapIt->second;
    double trackid = mapIt->first;
    if (energy > maxenergy){
      maxenergy = energy;
      objectTrack = trackid;
    }
  }

  return objectTrack;
}



int FDSelectionUtils::TrueParticleIDFromTotalRecoCharge(const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  // Make a map of the tracks which are associated with this object and the charge each contributes
  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int trackID = TrueParticleID(hit, rollup_unsaved_ids);
    trackMap[trackID] += hit->Integral();
  }

  // Pick the track with the highest charge as the 'true track'
  double highestCharge = 0;
  int objectTrack = -99999;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      objectTrack  = trackIt->first;
    }
  }
  return objectTrack;
}



int FDSelectionUtils::TrueParticleIDFromTotalRecoHits(const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  // Make a map of the tracks which are associated with this object and the number of hits they are the primary contributor to
  std::map<int,int> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    int trackID = TrueParticleID(hit, rollup_unsaved_ids);
    trackMap[trackID]++;
  }

  // Pick the track which is the primary contributor to the most hits as the 'true track'
  int objectTrack = -99999;
  int highestCount = -1;
  int NHighestCounts = 0;
  for (std::map<int,int>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
    if (trackIt->second > highestCount) {
      highestCount = trackIt->second;
      objectTrack  = trackIt->first;
      NHighestCounts = 1;
    }
    else if (trackIt->second == highestCount){
      NHighestCounts++;
    }
  }
  if (NHighestCounts > 1){
    std::cout<<"FDSelectionUtils::TrueParticleIDFromTotalRecoHits - There are " << NHighestCounts << " particles which tie for highest number of contributing hits (" << highestCount<<" hits).  Using FDSelectionUtils::TrueParticleIDFromTotalTrueEnergy instead."<<std::endl;
    objectTrack = FDSelectionUtils::TrueParticleIDFromTotalTrueEnergy(hits,rollup_unsaved_ids);
  }
  return objectTrack;
}

*/

double FDSelectionUtils::CompletenessFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& selected_hits, const std::vector<art::Ptr<recob::Hit> >& all_hits, int track_id){
  int num_matches_in_sel_hits = 0;
  int num_matches_in_all_hits = 0;
  for (unsigned int i_hit = 0; i_hit < selected_hits.size(); i_hit++){
    int matched_id = TruthMatchUtils::TrueParticleID(clockData, selected_hits[i_hit], 1); 
    if (matched_id==track_id) num_matches_in_sel_hits++;
  }
  for (unsigned int i_hit = 0; i_hit < all_hits.size(); i_hit++){
    int matched_id = TruthMatchUtils::TrueParticleID(clockData, all_hits[i_hit], 1); 
    if (matched_id==track_id) num_matches_in_all_hits++;
  }

  double completeness = 0;
  if (num_matches_in_all_hits > 0) completeness = 1.*num_matches_in_sel_hits/num_matches_in_all_hits;
  return completeness;
}

double FDSelectionUtils::HitPurityFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit> >& selected_hits, int track_id){
  int num_matches_in_sel_hits = 0;
  for (unsigned int i_hit = 0; i_hit < selected_hits.size(); i_hit++){
    int matched_id = TruthMatchUtils::TrueParticleID(clockData, selected_hits[i_hit], 1); 
    if (matched_id==track_id) num_matches_in_sel_hits++;
  }
  double hit_purity = 0;
  if (selected_hits.size() > 0) hit_purity = 1.*num_matches_in_sel_hits/selected_hits.size();
  return hit_purity;
}

bool FDSelectionUtils::IsInsideTPC(TVector3 position, double distance_buffer){
  double vtx[3] = {position.X(), position.Y(), position.Z()};
  bool inside = false;
  art::ServiceHandle<geo::Geometry> geom;
  geo::TPCID idtpc = geom->FindTPCAtPosition(vtx);

  if (geom->HasTPC(idtpc))
  {
    const geo::TPCGeo& tpcgeo = geom->GetElement(idtpc);
    double minx = tpcgeo.MinX(); double maxx = tpcgeo.MaxX();
    double miny = tpcgeo.MinY(); double maxy = tpcgeo.MaxY();
    double minz = tpcgeo.MinZ(); double maxz = tpcgeo.MaxZ();

    for (size_t c = 0; c < geom->Ncryostats(); c++)
    {
      const geo::CryostatGeo& cryostat = geom->Cryostat(c);
      for (size_t t = 0; t < cryostat.NTPC(); t++)
      {
        const geo::TPCGeo& tpcg = cryostat.TPC(t);
        if (tpcg.MinX() < minx) minx = tpcg.MinX();
        if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
        if (tpcg.MinY() < miny) miny = tpcg.MinY();
        if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
        if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
        if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
      }
    }

    //x
    double dista = fabs(minx - position.X());
    double distb = fabs(position.X() - maxx);
    if ((position.X() > minx) && (position.X() < maxx) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    //y
    dista = fabs(maxy - position.Y());
    distb = fabs(position.Y() - miny);
    if (inside && (position.Y() > miny) && (position.Y() < maxy) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    else inside = false;
    //z
    dista = fabs(maxz - position.Z());
    distb = fabs(position.Z() - minz);
    if (inside && (position.Z() > minz) && (position.Z() < maxz) &&
        (dista > distance_buffer) && (distb > distance_buffer)) inside = true;
    else inside = false;
  }

  return inside;
}

double FDSelectionUtils::CalculateTrackLength(const art::Ptr<recob::Track> track){
  double length = 0;
  if (track->NumberTrajectoryPoints()==1) return length; //Nothing to calculate if there is only one point

  for (size_t i_tp = 0; i_tp < track->NumberTrajectoryPoints()-1; i_tp++){ //Loop from the first to 2nd to last point
    TVector3 this_point(track->TrajectoryPoint(i_tp).position.X(),track->TrajectoryPoint(i_tp).position.Y(),track->TrajectoryPoint(i_tp).position.Z());
    if (!FDSelectionUtils::IsInsideTPC(this_point,0)){
      std::cout<<"FDSelectionUtils::CalculateTrackLength - Current trajectory point not in the TPC volume.  Skip over this point in the track length calculation"<<std::endl;
      continue;
    }
    TVector3 next_point(track->TrajectoryPoint(i_tp+1).position.X(),track->TrajectoryPoint(i_tp+1).position.Y(),track->TrajectoryPoint(i_tp+1).position.Z());
    if (!FDSelectionUtils::IsInsideTPC(next_point,0)){
      std::cout<<"FDSelectionUtils::CalculateTrackLength - Next trajectory point not in the TPC volume.  Skip over this point in the track length calculation"<<std::endl;
      continue;
    }

    length+=(next_point-this_point).Mag();
  }
  return length;
}



