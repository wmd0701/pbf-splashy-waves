#include <igl/writeOFF.h>
#include <thread>
#include "FluidSim.h"
#include "Gui.h"
#include "Simulator.h"

/*
 * This class is a GUI for our dummy simulation. It extends the basic GUI
 * defined in Gui.h. We could add more controls and visuals here, but we don't
 * need any additional functionality for this dummy simulation.
 */
class FluidGui : public Gui {
public:
	FluidSim *p_fluidSim = NULL;  // pointer to the dummy simulation

	FluidGui() {
		// create a new dummy simulation
		p_fluidSim = new FluidSim();

		// set this simulation as the simulation that is running in our GUI
		setSimulation(p_fluidSim);

		// start the GUI
		start();
	}

	virtual void updateSimulationParameters() override {
		// We don't have any simulation parameters to update periodically so we
		// don't need to do anything here
	};
};

int main(int argc, char *argv[]) {
	// create a new instance of the GUI for the dummy simulation
	new FluidGui();

	return 0;
}
