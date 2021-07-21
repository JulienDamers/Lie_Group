#ifndef SIVIA_H
#define SIVIA_H

#include "pyibex_SepProj.h"
#include "codac.h"
#include "codac-rob.h"
#include "ipegenerator.h"

void sivia(ibex::IntervalVector& map, pyibex::SepProj& fullSep, double epsilon);

void sivia_article(ibex::IntervalVector& map, pyibex::SepProj& fullSep, double epsilon, ipegenerator::Figure& fig);

#endif