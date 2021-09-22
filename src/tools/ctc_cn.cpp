//
// Created by julien-damers on 12/28/20.
//
#include "ctc_cn.h"

using namespace ibex;
using namespace codac;
using namespace std;

ctc_cn::ctc_cn(ContractorNetwork *cn, IntervalVectorVar *box)
        : Ctc(box->size()), m_cn(cn), m_box(box)
{

}

void ctc_cn::contract(ibex::IntervalVector &x)
{
    m_cn->reset_interm_vars();
    m_cn->contract({
                           {*m_box, x} // the "box" var is associated with "x"
                           // for this contraction procedure
                   });
}

void ctc_cn::contract(codac::TubeVector &x)
{
    int slice_id = 0;
    IntervalVector box_to_contract(x.size()+1);
    while (slice_id<x.nb_slices())
    {
        box_to_contract.put(0,x(slice_id));
        box_to_contract[box_to_contract.size()-1] = x[0].slice_tdomain(slice_id);
        this->contract(box_to_contract);
        if ( check_empty(box_to_contract))
        {
            x.set_empty();
            cout << "warning empty tube" << endl;
            break;
        }
        else
        {
            x.set(box_to_contract.subvector(0,box_to_contract.size()-2),slice_id);
        }
        slice_id++;
    }
}

bool ctc_cn::check_empty(ibex::IntervalVector &x)
{
    bool empty = false;
    int dim = 0;
    while(!empty && (dim<x.size()))
    {
        empty = x[dim].is_empty();
        dim++;
    }
    return(empty);
}