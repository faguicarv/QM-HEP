#include <fstream> // registrar datos
#include "Pythia8/Pythia.h"

int main() {

    int nevents = 1e4;

    Pythia8::Pythia pythia;

    pythia.readString("Beams:idA = -11"); // Beam de electrones
    pythia.readString("Beams:idB = 11"); // Beam de positrones
    pythia.readString("Beams:eCM = 91.1876"); // Z-Boson resonance energy

    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off"); // Apagamos canales de desintegración del Z
    pythia.readString("23:onIfMatch = 15 -15"); // Activamos sólo en canal de par taus

    pythia.init();

    std::ofstream output_file("e+e-_data.csv");
    output_file << "event,id,mass,px,py,pz,energy" << std::endl; //Encabezado archivo datos

    for(int i = 0; i < nevents; i++)
    {
        if(!pythia.next()) continue;

        for(int j = 0; j < pythia.event.size(); j++)
        {
            int id = pythia.event[j].id(); // variable guarda id partícula

            if(id == 23 || id == 15 || id == -15) // guardamos sólo partículas de interés (se puede definir con un bool también) ESTA FUNCIONA PERO PROBAMOS OTRA
            {

                double px = pythia.event[j].px(); // momentum partícula
                double py = pythia.event[j].py();
                double pz = pythia.event[j].pz();
                double mass = pythia.event[j].m();
                double E = pythia.event[j].e();

                output_file << i << ", "
                            << id << ", "
                            << mass << ", "
                            << px << ", "
                            << py << ", "
                            << pz << ", "
                            << E << std::endl;

                std::cout << id << " " << px << " " << py << " " << pz << std::endl;
            }
        }
    }

    output_file.close();
    return 0;
}
