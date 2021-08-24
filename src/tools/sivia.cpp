//
// Created by julien-damers on 1/15/21.
//
#include "sivia.h"
#include "figure.h"

using namespace ibex;
using namespace codac;
using namespace vibes;
using namespace pyibex;



void sivia(IntervalVector& map, Sep& Sep, double epsilon)
{

    int bisections = 0;

    stack<IntervalVector> s;
    s.push(map);


    while(!s.empty())
    {
        IntervalVector box = s.top();
        //DEBUG_MSG("X: " << box );
        s.pop();
        IntervalVector boxIn(box);
        IntervalVector boxOut(box);
        Sep.separate(boxIn,boxOut);
        //DEBUG_MSG("X_in: " << boxIn );
        //DEBUG_MSG("X_out: " << boxOut );
        //DEBUG_MSG(" ");
        if (boxOut[0].is_empty() && boxOut[1].is_empty())
        {
            drawBox(box[0].lb(),box[0].ub(), box[1].lb(),box[1].ub(), "#009E73[#56B4E9]");
        }
        else if (boxIn[0].is_empty() && boxIn[1].is_empty())
        {
            drawBox(box[0].lb(),box[0].ub(), box[1].lb(),box[1].ub(), "#D55E00[#CC79A7]");
        }
        else
        {
            if (box.max_diam()>epsilon)
            {
                int i = box.subvector(0,1).extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = box.bisect(i);
                s.push(p.first);
                s.push(p.second);
                bisections++;
            }
            else
            {
                drawBox(box[0].lb(),box[0].ub(),box[1].lb(),box[1].ub(),"#E69F00[#F0E442]");

            }
        }

    }
    cout << bisections << " bisectionss have been done" << endl;
}

void sivia(IntervalVector& map, Ctc& Ctc, double epsilon)
{

    int bisections = 0;

    stack<IntervalVector> s;
    s.push(map);


    while(!s.empty())
    {
        IntervalVector box = s.top();
        s.pop();
        IntervalVector boxOut(box);
        Ctc.contract(boxOut);
        if (boxOut[0].is_empty() && boxOut[1].is_empty())
        {
            drawBox(box[0].lb(),box[0].ub(), box[1].lb(),box[1].ub(), "b[#648FFF]");
        }
        else
        {
            if (box.max_diam()>epsilon)
            {
                int i = box.subvector(0,1).extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = box.bisect(i);
                s.push(p.first);
                s.push(p.second);
                bisections++;
            }
            else
            {
                drawBox(box[0].lb(),box[0].ub(),box[1].lb(),box[1].ub(),"y[#FFB000]");

            }
        }

    }
    cout << bisections << " bisections have been done" << endl;
}


void sivia_article(ibex::IntervalVector& map, ibex::Sep& Sep, double epsilon, ipegenerator::Figure& fig)
{
    stack<IntervalVector> s;
    s.push(map);

    int bisections = 0;

    while(!s.empty())
    {
        IntervalVector box = s.top();
        s.pop();
        IntervalVector boxIn(box);
        IntervalVector boxOut(box);
        Sep.separate(boxIn,boxOut);
        if (boxOut[0].is_empty() && boxOut[1].is_empty())
        {
            fig.draw_box(box,"colorBlindOutStroke","colorBlindOutFill");
        }
        else if (boxIn[0].is_empty() && boxIn[1].is_empty())
        {
            fig.draw_box(box,"colorBlindInStroke","colorBlindInFill");
        }
        else
        {
            if (box.max_diam()>epsilon)
            {
                int i = box.subvector(0,1).extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = box.bisect(i);
                s.push(p.first);
                s.push(p.second);
                bisections++;
            }
            else
            {
                fig.draw_box(box,"colorBlindMaybeStroke","colorBlindMaybeFill");

            }
        }

    }
    cout << bisections << " bisections have been done" << endl;

}

void sivia_article(ibex::IntervalVector& map, ibex::Ctc& Ctc, double epsilon, ipegenerator::Figure& fig)
{
    stack<IntervalVector> s;
    s.push(map);

    int bisections = 0;

    while(!s.empty())
    {
        IntervalVector box = s.top();
        s.pop();
        IntervalVector boxOut(box);
        Ctc.contract(boxOut);
        if (boxOut[0].is_empty() && boxOut[1].is_empty())
        {
            fig.draw_box(box,"blue","colorBlindOut");
        }
        else
        {
            if (box.max_diam()>epsilon)
            {
                int i = box.subvector(0,1).extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = box.bisect(i);
                s.push(p.first);
                s.push(p.second);
                bisections++;
            }
            else
            {
                fig.draw_box(box,"gold","colorBlindMaybe");

            }
        }

    }
    cout << bisections << " bisectionss have been done" << endl;

}

