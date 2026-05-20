#include <fstream> // registrar datos
#include "Pythia8/Pythia.h"

int main() {

    int nevents = 1e4;

    Pythia8::Pythia pythia;

    pythia.readString("Beams:idA = -11"); // Beam de electrones
    pythia.readString("Beams:idB = 11"); // Beam de positrones
    pythia.readString("Beams:eCM = 91.1876"); // Z-Boson resonance energy

    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    // pythia.readString("23:onMode = off"); // Apagamos canales de desintegración del Z
    pythia.readString("23:onMode = on"); // abrimos todos los canales de desintegración del Z
    // pythia.readString("23:onIfMatch = 15 -15"); // Activamos sólo en canal de par taus

    pythia.init();

    std::ofstream output_file("e+e-_data.csv"); // define output
    output_file << "event,id,mass,px,py,pz,energy" << std::endl; // Encabezado archivo datos

    for(int i = 0; i < nevents; i++) // Corremos de 0 a nevents-1
    {
        if(!pythia.next()) continue;

        for(int j = 0; j < pythia.event.size(); j++) // Vemos dentro de todos los eventos generados
        {
            int id = pythia.event[j].id(); // variable guarda id partícula
            int pos_daug1 = pythia.event[j].daughter1(); // Da la posición en el evento de la hija, no la id
            int pos_daug2 = pythia.event[j].daughter2();

            int id_daug1 = pythia.event[pos_daug1].id();
            int id_daug2 = pythia.event[pos_daug2].id();

            int status = pythia.event[j].status();


            if((id == 23) && std::abs(id_daug1) == 15)
            {
                output_file << i << ", " << id << ", " << pythia.event[j].m() << ", " << pythia.event[j].px() << ", " << pythia.event[j].py() << ", " << pythia.event[j].pz() << ", " << pythia.event[j].e() << std::endl; // Guardamos datos del Z madre

                output_file << i << ", " << id_daug1 << ", " << pythia.event[pos_daug1].m() << ", " << pythia.event[pos_daug1].px() << ", " << pythia.event[pos_daug1].py() << ", " << pythia.event[pos_daug1].pz() << ", " << pythia.event[pos_daug1].e() << std::endl; // Guardamos datos del tau 1


                output_file << i << ", " << id_daug2 << ", " << pythia.event[pos_daug2].m() << ", " << pythia.event[pos_daug2].px() << ", " << pythia.event[pos_daug2].py() << ", " << pythia.event[pos_daug2].pz() << ", " << pythia.event[pos_daug2].e() << std::endl; // Guardamos datos del tau 2
            }
        }

    }

    output_file.close();
    return 0;
}
