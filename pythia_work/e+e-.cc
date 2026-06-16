#include <fstream> // registrar datos
#include "Pythia8/Pythia.h"

int main(int argc, char* argv[]) {

    if (argc < 2)
    {
        std::cerr << "Uso: " << argv[0] << " <random_seed>" << std::endl;
        return 1;
    }

    int seed = std::stoi(argv[1]); // argumento a entero

    int nevents = 3e6;

    Pythia8::Pythia pythia;

    // Paralelizar pythia
    pythia.readString("Random:setSeed = on");  // Control manual de semilla
    pythia.readString("Random:seed = " + std::to_string(seed)); // Asigna semilla única

    pythia.readString("Beams:idA = -11"); // Beam de positrones
    pythia.readString("Beams:idB = 11"); // Beam de electrones
    pythia.readString("Beams:eCM = 91.1876"); // Z-Boson resonance energy

    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    // pythia.readString("23:onMode = on"); // abrimos todos los canales de desintegración del Z
    pythia.readString("23:onMode = off"); // Apagamos canales de desintegración del Z
    pythia.readString("23:onIfMatch = 15 -15"); // Activamos sólo el canal de par taus

    pythia.init();

//     std::ofstream output_file("e+e-_data.csv"); // define output SOLO PARA CORRER EN UN NODO
    // Output paralelizado
    std::string filename = "Z_decay_" + std::to_string(seed) + ".csv";
    std::ofstream output_file(filename);
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
                    // EN ESTE PUNTO HAY QUE LEER TODAS LAS NIETAS Y SELECCIONAR SÓLO LAS QUE NOS SIRVA
                    // EN VEZ DE BUSCAR LA ID UNA A UNA, MEJOR TENER LA LISTA COMPLETA Y DE AHÍ SELECCIONAR LO INTERESANTE

                    // Obtenemos los índices de las hijas del Tau 1
                    std::vector<int> grdaug_tau1 = pythia.event[pos_daug1].daughterList();
                    // Veamos ahora el Tau 2
                    std::vector<int> grdaug_tau2 = pythia.event[pos_daug2].daughterList();

                    // Para resolver los casos cuando se incluye un fotón radiado o estados intermedios que emiten los mismos productos deseados (como el rho) vamos a buscar eventos que sólo tengan dos eventos, porque es importante mantener sólo un canal de interés en este trabajo
                    if(grdaug_tau1.size() == 2 && grdaug_tau2.size() == 2)
                    {

                        // Recorremos el vector en busca de las ID deseadas
                        int idx_pi1 = -1; // Variables auxiliares para buscar IDs
                        int idx_nu1 = -1;

                        for(int idx : grdaug_tau1) // idx son las posiciones donde hay hijas. Hay que calzar idx con piones y neutrinos
                        {
                            int aux = std::abs(pythia.event[idx].id());
                            if(aux == 211) // Si es pion
                            {
                                idx_pi1 = idx; // Guardamos el pion

                            }
                            else if(aux == 16) // Si es neutrino
                            {
                                idx_nu1 = idx; // Guardamos el neutrino
                            }
                        } // Hasta acá recorrimos el Tau 1 y guardamos las variables.


                        // Recorremos el vector en busca de las ID deseadas
                        int idx_pi2 = -1; // Variables auxiliares para buscar IDs
                        int idx_nu2 = -1;
                        for(int idx : grdaug_tau2) // idx son las posiciones donde hay hijas. Hay que calzar idx con piones y neutrinos
                        {
                            int aux = std::abs(pythia.event[idx].id());
                            if(aux == 211) // Si es pion
                            {
                                idx_pi2 = idx; // Guardamos el pion
                            }
                            else if(aux == 16) // Si es neutrino
                            {
                                idx_nu2 = idx; // Guardamos el neutrino
                            }
                        }

                        if(idx_pi1 != -1 && idx_nu1 != -1 && idx_pi2 != -1 && idx_nu2 != -1)
                        {
                            // Si pasaron todos estos condicionales entonces encontramos el canal objetivo y que lo diga en la consola
                            std::cout << "¡Canal Z-boson -> Tau/Antitau -> pi^ nu_tau detectado en evento " << i << "!" << std::endl;

                            // Guardamos información en el output
                            output_file << i << ", " << id << ", " << pythia.event[j].m() << ", " << pythia.event[j].e() << ", " << pythia.event[j].px() << ", " << pythia.event[j].py() << ", " << pythia.event[j].pz() << std::endl; // Guardamos datos del Z madre

                            output_file << i << ", " << id_daug1 << ", " << pythia.event[pos_daug1].m() << ", " << pythia.event[pos_daug1].e() << ", " << pythia.event[pos_daug1].px() << ", " << pythia.event[pos_daug1].py() << ", " << pythia.event[pos_daug1].pz() << std::endl; // Guardamos datos del Tau 1

                            output_file << i << ", " << id_daug2 << ", " << pythia.event[pos_daug2].m() << ", " << pythia.event[pos_daug2].e() << ", " << pythia.event[pos_daug2].px() << ", " << pythia.event[pos_daug2].py() << ", " << pythia.event[pos_daug2].pz() << std::endl; // Guardamos datos del Tau 2

                            output_file << i << ", " << pythia.event[idx_pi1].id() << ", " << pythia.event[idx_pi1].m() << ", " <<  pythia.event[idx_pi1].e() << ", " << pythia.event[idx_pi1].px() << ", " << pythia.event[idx_pi1].py() << ", " << pythia.event[idx_pi1].pz() << std::endl; // Guardamos nieta 1 del Tau 1

                            output_file << i << ", " << pythia.event[idx_nu1].id() << ", " << pythia.event[idx_nu1].m() << ", " <<  pythia.event[idx_nu1].e() << ", " << pythia.event[idx_nu1].px() << ", " << pythia.event[idx_nu1].py() << ", " << pythia.event[idx_nu1].pz() << std::endl; // Guardamos nieta 2 del Tau 1

                            output_file << i << ", " << pythia.event[idx_pi2].id() << ", " << pythia.event[idx_pi2].m() << ", " <<  pythia.event[idx_pi2].e() << ", " << pythia.event[idx_pi2].px() << ", " << pythia.event[idx_pi2].py() << ", " << pythia.event[idx_pi2].pz() << std::endl; // Guardamos nieta 1 del Tau 2

                            output_file << i << ", " << pythia.event[idx_nu2].id() << ", " << pythia.event[idx_nu2].m() << ", " <<  pythia.event[idx_nu2].e() << ", " << pythia.event[idx_nu2].px() << ", " << pythia.event[idx_nu2].py() << ", " << pythia.event[idx_nu2].pz() << std::endl; // Guardamos nieta 2 del Tau 2

                        }
                    }
                }
            }
        }
    }

    output_file.close();
    return 0;
}
