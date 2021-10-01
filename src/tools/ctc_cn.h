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
            ctc_cn(codac::ContractorNetwork *cn, codac::IntervalVectorVar *box);
            void contract(ibex::IntervalVector& x);
            //void contract(codac::TubeVector& x);
            bool check_empty(ibex::IntervalVector& x);

        private:
            codac::ContractorNetwork* m_cn;
            codac::IntervalVectorVar* m_box;
    };

#endif //PAPER_LIE_GROUP_CTC_CN_H
