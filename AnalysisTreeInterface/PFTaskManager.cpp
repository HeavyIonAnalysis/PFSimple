#include "PFTaskManager.h"

#include "SimpleFinder.h"

void PFTaskManager::Run(long long int nEvents) {
  std::cout << "PFTaskManager::Run" << std::endl;
  nEvents = nEvents < 0 || nEvents > chain_->GetEntries() ? chain_->GetEntries() : nEvents;

  for (long long iEvent = 0; iEvent < nEvents; ++iEvent) {
    chain_->GetEntry(iEvent);

    if (iEvent % 100 == 0)
      std::cout << "Event # " << iEvent << " out of " << nEvents << "\r" << std::flush;

    auto* converter_in = (ConverterIn*) tasks_.at(kInConverter);
    SetDecay(converter_in->GetDecay());
    SetCuts(converter_in->GetCuts());

    SimpleFinder FCFinder;
    FCFinder.Init(converter_in->CreateInputContainer());
    FCFinder.SetDecay(converter_in->GetDecay());
    FCFinder.SetCuts(converter_in->GetCuts());
    FCFinder.SortTracks();
    FCFinder.FindParticles();
    auto* converter_out = ((ConverterOut*) tasks_.at(kOutConverter));
    converter_out->SetDecay(converter_in->GetDecay());
    converter_out->SetCandidates(FCFinder.GetMotherCandidates());
    converter_out->Exec();
    //     out_tree_->Fill();
  }// Event loop
}
