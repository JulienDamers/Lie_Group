//
// Created by julien-damers on 1/15/21.
//

#include "sivia.h"
//#include "figure.h"

using namespace ibex;
using namespace codac;
using namespace vibes;
using namespace pyibex;



void sivia(IntervalVector& map, SepProj& fullSep, double epsilon)
{

    int boxes_counter = 1;
    int bissection = 0;

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
            drawBox(box[0].lb(),box[0].ub(), box[1].lb(),box[1].ub(), "b[#648FFF]");
        }
        else if (boxIn[0].is_empty() && boxIn[1].is_empty())
        {
            drawBox(box[0].lb(),box[0].ub(), box[1].lb(),box[1].ub(), "r[#DC267F]");
        }
        else
        {
            if (box.max_diam()>epsilon)
            {
                int i = box.subvector(0,1).extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = box.bisect(i);
                s.push(p.first);
                s.push(p.second);
                bissection++;
            }
            else
            {
                drawBox(box[0].lb(),box[0].ub(),box[1].lb(),box[1].ub(),"y[#FFB000]");

            }
        }

    }
    cout << bissection << " bissections have been done" << endl;
}


void sivia_article(ibex::IntervalVector& map, pyibex::SepProj& fullSep, double epsilon, ipegenerator::Figure& fig)
{
    stack<IntervalVector> s;
    s.push(map);

    int bissection = 0;

    while(!s.empty())
    {
        IntervalVector box = s.top();
        s.pop();
        IntervalVector boxIn(box);
        IntervalVector boxOut(box);
        fullSep.separate(boxIn,boxOut);
        if (boxOut[0].is_empty() && boxOut[1].is_empty())
        {
            fig.draw_box(box,"blue","colorBlindOut");
        }
        else if (boxIn[0].is_empty() && boxIn[1].is_empty())
        {
            fig.draw_box(box,"red","colorBlindIn");
        }
        else
        {
            if (box.max_diam()>epsilon)
            {
                int i = box.subvector(0,1).extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = box.bisect(i);
                s.push(p.first);
                s.push(p.second);
                bissection++;
            }
            else
            {
                fig.draw_box(box,"gold","colorBlindMaybe");

            }
        }

    }
    cout << bissection << " bissections have been done" << endl;

}
