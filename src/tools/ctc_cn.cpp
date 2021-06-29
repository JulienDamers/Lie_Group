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
    for (int i = 0; i < m_intermediary_iv->size(); i++)
    {
        //cout << "intermediary vector[" << i << "] = " << *(*m_intermediary_iv)[i]<< endl;
        (*m_intermediary_iv)[i]->init(Interval(Interval::ALL_REALS));
    }
    for (int i = 0; i <  m_intermediary_i->size(); i++)
    {
        *(*m_intermediary_i)[i] = Interval(Interval::ALL_REALS);
    }
    m_cn->trigger_all_contractors();
    m_cn->contract();
    if (m_box_to_contract_cn->is_empty())
    {
        x.init(Interval(Interval::EMPTY_SET));
    }
    else
    {
        x.put(0,*m_box_to_contract_cn);
    }

}
