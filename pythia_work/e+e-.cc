#include <fstream> // registrar datos
#include "Pythia8/Pythia.h"

int main() {

    int nevents = 5e6;

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

            if(id == 23) // Sólo haremos el proceso para bosones Z
            {

                // Interrogamos a la primera hija del Z
                int pos_daug1 = pythia.event[j].daughter1(); // Da la posición en el evento de la hija, no la id
                int id_daug1 = pythia.event[pos_daug1].id(); // ID hija 1

                // Interrogamos a la segunda hija del Z
                int pos_daug2 = pythia.event[j].daughter2(); // Posición hija 2
                int id_daug2 = pythia.event[pos_daug2].id(); // ID hija 2

                if(std::abs(id_daug1) == 15 && std::abs(id_daug2) == 15) // Seguimos pero sólo si las hijas son un tau o antitau
                {
                    // Interrogamos al tau en busca de pi^+ y nu_tau
                    int pos_grdaug1_1 = pythia.event[pos_daug1].daughter1();
                    int pos_grdaug1_2 = pythia.event[pos_daug1].daughter2();

                    int id_grdaug1_1 = pythia.event[pos_grdaug1_1].id(); // Id nieta 1 del Z (producto del Tau)
                    int id_grdaug1_2 = pythia.event[pos_grdaug1_2].id(); // Id nieta 2 del Z (producto del Tau)

                    // Interrogamos al Tau 2 en busca de pi y nu_tau
                    int pos_grdaug2_1 = pythia.event[pos_daug2].daughter1();
                    int pos_grdaug2_2 = pythia.event[pos_daug2].daughter2();

                    int id_grdaug2_1 = pythia.event[pos_grdaug2_1].id(); // Id nieta 1 del Z (producto del Tau)
                    int id_grdaug2_2 = pythia.event[pos_grdaug2_2].id(); // Id nieta 2 del Z (producto del Tau)

                    if((std::abs(id_grdaug1_1) == 211 || std::abs(id_grdaug1_1) == 16) && (std::abs(id_grdaug1_2) == 211 || std::abs(id_grdaug1_2) == 16))
                    // Debemos buscar que id_grdaug1,2 sean pi (id=+/-211) y neutrino tau (id=+/-16)
                    {
                        if((std::abs(id_grdaug2_1) == 211 || std::abs(id_grdaug2_1) == 16) && (std::abs(id_grdaug2_2) == 211 || std::abs(id_grdaug2_2) == 16))
                        {


                            // Si pasaron todos estos condicionales entonces encontramos el canal objetivo y que lo diga en la consola
                            std::cout << "¡Canal Z-boson -> Tau/Antitau -> pi^ nu_tau detectado en evento " << i << "!" << std::endl;

                            output_file << i << ", " << id << ", " << pythia.event[j].m() << ", " << pythia.event[j].e() << ", " << pythia.event[j].px() << ", " << pythia.event[j].py() << ", " << pythia.event[j].pz() << std::endl; // Guardamos datos del Z madre

                            output_file << i << ", " << id_daug1 << ", " << pythia.event[pos_daug1].m() << ", " << pythia.event[pos_daug1].e() << ", " << pythia.event[pos_daug1].px() << ", " << pythia.event[pos_daug1].py() << ", " << pythia.event[pos_daug1].pz() << std::endl; // Guardamos datos del Tau 1

                            output_file << i << ", " << id_daug2 << ", " << pythia.event[pos_daug2].m() << ", " << pythia.event[pos_daug2].e() << ", " << pythia.event[pos_daug2].px() << ", " << pythia.event[pos_daug2].py() << ", " << pythia.event[pos_daug2].pz() << std::endl; // Guardamos datos del Tau 2

                            output_file << i << ", " << id_grdaug1_1 << ", " << pythia.event[pos_grdaug1_1].m() << ", " <<  pythia.event[pos_grdaug1_1].e() << ", " << pythia.event[pos_grdaug1_1].px() << ", " << pythia.event[pos_grdaug1_1].py() << ", " << pythia.event[pos_grdaug1_1].pz() << std::endl; // Guardamos nieta 1 del Tau 1

                            output_file << i << ", " << id_grdaug2_1 << ", " << pythia.event[pos_grdaug2_1].m() << ", " <<  pythia.event[pos_grdaug2_1].e() << ", " << pythia.event[pos_grdaug2_1].px() << ", " << pythia.event[pos_grdaug2_1].py() << ", " << pythia.event[pos_grdaug2_1].pz() << std::endl; // Guardamos nieta 1 del Tau 2

                            output_file << i << ", " << id_grdaug1_2 << ", " << pythia.event[pos_grdaug1_2].m() << ", " << pythia.event[pos_grdaug1_2].e() << ", " << pythia.event[pos_grdaug1_2].px() << ", " << pythia.event[pos_grdaug1_2].py() << ", " << pythia.event[pos_grdaug1_2].pz() << std::endl; // Guardamos nieta 2 del Tau 1

                            output_file << i << ", " << id_grdaug2_2 << ", " << pythia.event[pos_grdaug2_2].m() << ", " << pythia.event[pos_grdaug2_2].e() << ", " << pythia.event[pos_grdaug2_2].px() << ", " << pythia.event[pos_grdaug2_2].py() << ", " << pythia.event[pos_grdaug2_2].pz() << std::endl; // Guardamos nieta 2 del Tau 2

                        }
                    }
                }


            }
        }
    }

    output_file.close();
    return 0;
}
