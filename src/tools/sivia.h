#ifndef SIVIA_H
#define SIVIA_H

#include "pyibex_SepProj.h"
#include "codac.h"
#include "codac-rob.h"
#include "ipegenerator.h"

void sivia(ibex::IntervalVector& map, ibex::Sep& Sep, double epsilon);
void sivia(ibex::IntervalVector& map, ibex::Ctc& Ctc, double epsilon);

void sivia_article(ibex::IntervalVector& map, ibex::Sep& Sep, double epsilon, ipegenerator::Figure& fig);
void sivia_article(ibex::IntervalVector& map, ibex::Ctc& Ctc, double epsilon, ipegenerator::Figure& fig);
#endif