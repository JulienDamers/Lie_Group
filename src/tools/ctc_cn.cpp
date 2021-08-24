//
// Created by julien-damers on 12/28/20.
//
#include "ctc_cn.h"

using namespace ibex;
using namespace codac;
using namespace std;

ctc_cn::ctc_cn(ContractorNetwork* cn, IntervalVector* box_to_contract_cn,
               vector<IntervalVector*>* intermediary_iv,
               vector<Interval*>* intermediary_i): Ctc(box_to_contract_cn->size())
{
    m_cn = cn;
    m_box_to_contract_cn = box_to_contract_cn;
    m_intermediary_iv = intermediary_iv;
    m_intermediary_i = intermediary_i;
}

void ctc_cn::contract(ibex::IntervalVector &x)
{
    m_box_to_contract_cn->clear();
    m_box_to_contract_cn->put(0,x);
    for (size_t i = 0; i < m_intermediary_iv->size(); i++)
    {
        //cout << "intermediary interval vector[" << i << "] = " << *(*m_intermediary_iv)[i]<< endl;
        (*m_intermediary_iv)[i]->init(Interval::ALL_REALS);
    }
    for (size_t i = 0; i <  m_intermediary_i->size(); i++)
    {
        *(*m_intermediary_i)[i] = Interval::ALL_REALS;
        //cout << "intermediary interval[" << i << "] = " << *(*m_intermediary_i)[i]<< endl;
    }
    m_cn->trigger_all_contractors();
    m_cn->contract();
    if ( check_empty(*m_box_to_contract_cn))
    {
        x.init(Interval::EMPTY_SET);
    }
    else
    {
        x.put(0,*m_box_to_contract_cn);
    }

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