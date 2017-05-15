#ifndef PLANT_PROPAGATION_H
#define PLANT_PROPAGATION_H

#include "framework.h"

Population runPPA(Parameters *params, ValueCollector &vc);
Population runPPALevy(Parameters *params, ValueCollector &vc);

#endif