#include "Manager.h"

void Manager::Run(int n_events)
{
  InConverter converter_in;
  converter_in.SetTrackCuts(track_cuts_);
  converter_in.InitAnalysisTree(in_file_name_, in_tree_name_);
  converter_in.SetIsShine(is_shine_);

  OutConverter converter_out;
  converter_out.SetOutFileName("KFPS_output.root");
  converter_out.InitAT();

  std::cout << converter_in.GetNEvents() << std::endl;

  if(n_events<0 || n_events>converter_in.GetNEvents())
    n_events = converter_in.GetNEvents();
  
  for(int iEvent=0; iEvent<n_events; iEvent++)
  {
    InputContainer inputInfo = converter_in.CreateInputContainer(iEvent);
    
    SimpleFinder FCFinder;
    FCFinder.Init(inputInfo);
    FCFinder.SetCuts(cuts_);
    FCFinder.SortTracks();
    FCFinder.FindParticles();

    converter_out.WriteCandidates(FCFinder.GetLambdas());
  }
  
  converter_out.Finish();
}