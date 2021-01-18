//
// Created by julien-damers on 08/01/2020.
//

#ifndef LIEGROUPTOOLS_LIE_GROUP_EX4_SEPARATOR_H
#define LIEGROUPTOOLS_LIE_GROUP_EX4_SEPARATOR_H

#include "tubex.h"
#include "tubex-rob.h"
#include "ibex.h"

class lie_group_ex4_separator : public ibex::Sep
{
    public:
        lie_group_ex4_separator(tubex::TubeVector* a, ibex::IntervalVector* constraint);
        ~lie_group_ex4_separator();
        void separate(ibex::IntervalVector& Xin, ibex::IntervalVector& Xout);

    private:
        ibex::SepFwdBwd* sepPhi;
        ibex::IntervalVector* constraintVector;
        ibex::Function* fy;
        tubex::TubeVector* reference;
};


#endif //LIEGROUPTOOLS_LIE_GROUP_EX4_SEPARATOR_H
