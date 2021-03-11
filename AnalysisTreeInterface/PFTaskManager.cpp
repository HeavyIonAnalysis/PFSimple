#include "PFTaskManager.h"

#include "SimpleFinder.h"

void PFTaskManager::Run(long long int nEvents) {
  std::cout << "PFTaskManager::Run" << std::endl;

  auto* chain = man_->GetChain();

  nEvents = nEvents < 0 || nEvents > chain->GetEntries() ? chain->GetEntries() : nEvents;

  for (long long iEvent = 0; iEvent < nEvents; ++iEvent) {
    chain->GetEntry(iEvent);

    if (iEvent % 100 == 0)
      std::cout << "Event # " << iEvent << " out of " << nEvents << "\r" << std::flush;

    auto& tasks = man_->Tasks();

    auto* converter_in = (ConverterIn*) tasks.at(kInConverter);
    SetDecay(converter_in->GetDecay());
    SetCuts(converter_in->GetCuts());
//
    SimpleFinder FCFinder;
    FCFinder.Init(converter_in->CreateInputContainer());
    FCFinder.SetDecay(converter_in->GetDecay());
    FCFinder.SetCuts(converter_in->GetCuts());
    FCFinder.SortTracks();
    FCFinder.FindParticles();
//
    auto* converter_out = ((ConverterOut*) tasks.at(kOutConverter));
    converter_out->SetDecay(converter_in->GetDecay());
    converter_out->SetCandidates(FCFinder.GetMotherCandidates());
    converter_out->Exec();
    //     out_tree_->Fill(); b
  }// Event loop
}
