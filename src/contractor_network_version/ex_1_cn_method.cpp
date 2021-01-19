//
// Created by julien-damers on 04/01/2021.
//

#include "tubex.h"
#include "tubex-rob.h"
#include "tubex-capd.h"
#include "tools.h"
#include "chrono"
#include "ctc_cn.h"



using namespace std;
using namespace tubex;
using namespace ibex;
using namespace vibes;
using namespace pyibex;

int main(int argc, char* argv[])
{

    Tube::enable_syntheses();

    Interval domain(0,8);
    double timestep = 0.01;
    IntervalVector x0({{0,0},{1,1}});
    Function f("x","y","(1;-y)");

    // Create the reference curve using CAPD
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);


    double epsilon = 0.1;
    IntervalVector X0({{0,1},{2,3}});
    ibex::Function phi("x1", "x2", "t", "a1", "a2", "(x1 - a1; x2/a2)");
    ibex::CtcFwdBwd ctc_phi(phi,X0);
    ibex::CtcNotIn ctc_phi_not(phi,X0);
    tubex::CtcEval ctc_eval;


    // Complete version

    ContractorNetwork cn_out;
    vector<IntervalVector*> intermediary_iv_out;
    IntervalVector box_out = IntervalVector(3);
    IntervalVector z_out = IntervalVector(2); // a(t)
    intermediary_iv_out.push_back(&z_out);
    cn_out.add(ctc_eval,{box_out[2],z_out,a}); // t = x3
    cn_out.add(ctc_phi,{box_out,z_out});
    ctc_cn ctcCn_out(&cn_out,&box_out,&intermediary_iv_out);

    ContractorNetwork cn_in;
    vector<IntervalVector*> intermediary_iv_in;
    IntervalVector box_in = IntervalVector(3);
    IntervalVector z_in = IntervalVector(2); // a(t)
    intermediary_iv_in.push_back(&z_in);
    cn_in.add(ctc_eval,{box_in[2],z_in,a}); // t = x3
    cn_in.add(ctc_phi_not,{box_in,z_in});
    ctc_cn ctcCn_in(&cn_in,&box_in,&intermediary_iv_in);

    SepCtcPair sep(ctcCn_in,ctcCn_out);
    SepProj proj(sep,Interval(0,4),epsilon);

    IntervalVector map({{0,4},{-0.2,4}});

    beginDrawing();
    VIBesFigMap fig("ex_1_cn_method");
    fig.set_properties(50, 50, 800, 800);
    fig.axis_limits(map.subvector(0,1));
    stack<IntervalVector> s;


    auto start = chrono::steady_clock::now();

    // Separator Version

    s.push(map);
    while ( !s.empty())
    {
        IntervalVector current_box = s.top();
        s.pop();
        IntervalVector boxIn = current_box;
        IntervalVector boxOut = current_box;
        proj.separate(boxIn,boxOut);
        if (boxOut[0].is_empty() && boxOut[1].is_empty())
        {
            drawBox(current_box,"k[b]");
        }
        else if(boxIn[0].is_empty() && boxIn[1].is_empty())
        {
            drawBox(current_box, "k[r]");
        }
        else
        {
            if (current_box.max_diam()>epsilon)
            {
                int i = current_box.extr_diam_index(false);
                pair<IntervalVector, IntervalVector> p = current_box.bisect(i);
                s.push(p.first);
                s.push(p.second);
            }
            else
            {
                drawBox(current_box,"k[y]");
            }

        }

    }

    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" << endl;
    fig.add_tube(&a,"reference",0,1);
    drawBox(X0,"g[]");
    fig.axis_limits(map);
    fig.show();











}