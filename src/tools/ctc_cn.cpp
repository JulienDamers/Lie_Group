//
// Created by julien-damers on 12/28/20.
//
#include "ctc_cn.h"

using namespace ibex;
using namespace codac;
using namespace std;

CtcCn::CtcCn(ContractorNetwork *cn, IntervalVectorVar *box)
        : Ctc(box->size()), m_cn(cn), m_box(box)
{

}

void CtcCn::contract(ibex::IntervalVector &x)
{
    m_cn->reset_interm_vars();
    m_cn->contract({
                           {*m_box, x} // the "box" var is associated with "x"
                           // for this contraction procedure
                   });
}