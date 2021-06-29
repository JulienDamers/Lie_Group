//
// Created by julien-damers on 18/12/2019.
//



#ifndef LIEGROUPTOOLS_LIE_GROUP_EX3_SEPARATOR_H
#define LIEGROUPTOOLS_LIE_GROUP_EX3_SEPARATOR_H

#include "codac.h"
#include "codac-rob.h"
#include "ibex.h"


    class lie_group_ex3_separator : public ibex::Sep
{
    public:
        lie_group_ex3_separator(codac::TubeVector* a, ibex::IntervalVector* constraint);
        ~lie_group_ex3_separator();
        void separate(ibex::IntervalVector &Xin, ibex::IntervalVector &Xout);



    private:
        codac::TubeVector* reference;
        ibex::IntervalVector* constraintVector;
        ibex::Interval constraintDomain;
        ibex::Function* fy;
        ibex::Function* fdom;
        ibex::SepFwdBwd* sepPhi;
        ibex::SepFwdBwd* sepDom;
        ibex::Sep* fullSep;


    };


#endif //LIEGROUPTOOLS_LIE_GROUP_EX3_SEPARATOR_H
