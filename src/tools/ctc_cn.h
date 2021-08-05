//
// Created by julien-damers on 12/28/20.
//

#ifndef PAPER_LIE_GROUP_CTC_CN_H
#define PAPER_LIE_GROUP_CTC_CN_H

#include "ibex.h"
#include "codac.h"

    class ctc_cn : public ibex::Ctc
    {
        public:
            ctc_cn(codac::ContractorNetwork* cn, ibex::IntervalVector* box_to_contract_cn,
                   std::vector<ibex::IntervalVector*>* intermediary_iv = new std::vector<ibex::IntervalVector*>(),
                   std::vector<ibex::Interval*>* intermediary_i = new std::vector<ibex::Interval*>());
            void contract(ibex::IntervalVector& x);
            void contract(codac::TubeVector& x);
            bool check_empty(ibex::IntervalVector& x);

        private:
            codac::ContractorNetwork* m_cn;
            ibex::IntervalVector* m_box_to_contract_cn;
            std::vector<ibex::IntervalVector*>* m_intermediary_iv;
            std::vector<ibex::Interval*>* m_intermediary_i;
    };

#endif //PAPER_LIE_GROUP_CTC_CN_H
