#ifndef SIVIA_H
#define SIVIA_H

#include "pyibex_SepProj.h"
#include "tubex.h"
#include "tubex-rob.h"

void sivia(ibex::IntervalVector& map, pyibex::SepProj fullSep, double epsilon);

#endif