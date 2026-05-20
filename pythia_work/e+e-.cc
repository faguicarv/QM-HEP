#include <fstream> // registrar datos
#include "Pythia8/Pythia.h"

int main() {

    int nevents = 1e6;

    Pythia8::Pythia pythia;

    pythia.readString("Beams:idA = -11"); // Beam de positrones
    pythia.readString("Beams:idB = 11"); // Beam de electrones
    pythia.readString("Beams:eCM = 91.1876"); // Z-Boson resonance energy

    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    // pythia.readString("23:onMode = off"); // Apagamos canales de desintegración del Z
    pythia.readString("23:onMode = on"); // abrimos todos los canales de desintegración del Z
    // pythia.readString("23:onIfMatch = 15 -15"); // Activamos sólo el canal de par taus

    pythia.init();

    std::ofstream output_file("e+e-_data.csv"); // define output
    output_file << "event,id,mass,energy,px,py,pz" << std::endl; // Encabezado archivo datos

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

            int pos_grdaug1 = pythia.event[pos_daug1].daughter1(); // Busquemos eventos en que el tau desintegre en pi+ y neutrino tau
            int pos_grdaug2 = pythia.event[pos_daug1].daughter2();
            int id_grdaug1 = pythia.event[pos_grdaug1].id(); // Id del producto de desintegración del tau
            int id_grdaug2 = pythia.event[pos_grdaug2].id();

            // Debemos buscar que id_grdaug1,2 sean pi+ (id=211) y neutrino tau (id=16)

            int status = pythia.event[j].status();


            if((id == 23) && std::abs(id_daug1) == 15)
            {
                if((id_grdaug1 == 211 || id_grdaug1 == 16) && (id_grdaug2 == 211 || id_grdaug2 == 16))
                {

                    output_file << i << ", " << id << ", " << pythia.event[j].m() << ", " << pythia.event[j].e() << ", " << pythia.event[j].px() << ", " << pythia.event[j].py() << ", " << pythia.event[j].pz() << std::endl; // Guardamos datos del Z madre

                    output_file << i << ", " << id_daug1 << ", " << pythia.event[pos_daug1].m() << ", " << pythia.event[pos_daug1].e() << ", " << pythia.event[pos_daug1].px() << ", " << pythia.event[pos_daug1].py() << ", " << pythia.event[pos_daug1].pz() << std::endl; // Guardamos datos del tau 1


                    output_file << i << ", " << id_daug2 << ", " << pythia.event[pos_daug2].m() << ", " << pythia.event[pos_daug2].e() << ", " << pythia.event[pos_daug2].px() << ", " << pythia.event[pos_daug2].py() << ", " << pythia.event[pos_daug2].pz() << std::endl; // Guardamos datos del tau 2

                    output_file << i << ", " << id_grdaug1 << ", " << pythia.event[pos_grdaug1].m() << ", " <<  pythia.event[pos_grdaug1].e() << ", " << pythia.event[pos_grdaug1].px() << ", " << pythia.event[pos_grdaug1].py() << ", " << pythia.event[pos_grdaug1].pz() << std::endl; // Guardar datos tau o pi+

                    output_file << i << ", " << id_grdaug2 << ", " << pythia.event[pos_grdaug2].m() << ", " << pythia.event[pos_grdaug2].e() << ", " << pythia.event[pos_grdaug2].px() << ", " << pythia.event[pos_grdaug2].py() << ", " << pythia.event[pos_grdaug2].pz() << std::endl; // Guardar datos tau o pi+

                    std::cout << "¡Canal Z-boson -> Tau -> pi^+ nu detectado en evento " << i << "!" << std::endl;
                }
            }
        }

    }

    output_file.close();
    return 0;
}
