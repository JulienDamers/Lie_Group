//
// Created by julien-damers on 1/15/21.
//

#include "sivia.h"

using namespace ibex;
using namespace tubex;
using namespace vibes;
using namespace pyibex;


void sivia(IntervalVector& map, SepProj fullSep, double epsilon)
{
    stack<IntervalVector> s;
    s.push(map);

    while(!s.empty())
    {
        IntervalVector box = s.top();
        s.pop();
        IntervalVector boxIn(box);
        IntervalVector boxOut(box);
        fullSep.separate(boxIn,boxOut);
        if (boxOut[0].is_empty() && boxOut[1].is_empty())
        {
            drawBox(box[0].lb(),box[0].ub(), box[1].lb(),box[1].ub(), "b[c]");
        }
        else if (boxIn[0].is_empty() && boxIn[1].is_empty())
        {
            drawBox(box[0].lb(),box[0].ub(), box[1].lb(),box[1].ub(), "r[m]");
        }
        else
        {
            if (box.max_diam()>epsilon)
            {
                int i = box.subvector(0,1).extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = box.bisect(i);
                s.push(p.first);
                s.push(p.second);
            }
            else
            {
                drawBox(box[0].lb(),box[0].ub(),box[1].lb(),box[1].ub(),"w[w]");

            }
        }

    }
}