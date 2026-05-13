#include <iostream>

#include "Pythia8/Pythia.h"

int main()
{
	int nevents = 1e4;

	Pythia8::Pythia pythia;

	pythia.readString("Beams:idA = 2212");  // Beam primario (izq) y su id (proton)
	pythia.readString("Beams:idB = 2212");  // Beam primario (der) de proton
	pythia.readString("Beams:eCM = 14.e3"); // Energía del centro de masa (Pythia tiene por defecto GeV para energías)

	pythia.readString("HiggsSM:ffbar2H = on"); // Higgs por fusión de gluones
    pythia.readString("25:onMode = off");  // Apagar canales desintegración del higgs (id 25)
    pythia.readString("25:onIfMatch = 15 -15"); // Encender solamente tau/antitau (15/-15) provenientes de la desintegración del higgs
    pythia.readString("15:onIfMatch = 16 211"); // si tenemos un tau+, seleccionamos el canal de desintegración de neutrino tau (16) y pion + (211)
	//pythia.readString("SoftQCD:all = on");
	//pythia.readString("HardQCD:all = on");

	Pythia8::Hist hpz("Momentum Distribution", 100, -10, 10);

	pythia.init();

	for(int i = 0; i < nevents; i++)
	{
		if(!pythia.next()) continue;

		int entries = pythia.event.size();

		std::cout << "Event: "<< std::endl;
		std::cout << "Event size: "<< entries << std::endl;

        for(int j = 0; j < entries; j++)
		{
            int id = pythia.event[j].id();

			double m = pythia.event[j].m();

			double px = pythia.event[j].px();
			double py = pythia.event[j].py();
			double pz = pythia.event[j].pz();

			double pabs = sqrt(pow(px, 2) + pow(py,2) + pow(pz,2));

			hpz.fill(pz);

			std::cout << id << " " << m << " " << pabs << std::endl;
        }
    }

	std::cout << hpz << std::endl;

	Pythia8::HistPlot hpl("higgsdecay");
	hpl.frame("output", "Momentum Distribution", "Momentum", "Entries");
	hpl.add(hpz);
	hpl.plot();

	return 0;
}
