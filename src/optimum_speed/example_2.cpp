//
// Created by julien-damers on 12/12/2019.
//
#include <iostream>
#include <random>
#include "tools.h"
#include "codac.h"
#include "codac-capd.h"

#include "chrono"

using namespace std;
using namespace ibex;
using namespace codac;
using namespace vibes;
using namespace pyibex;





int main (int argc, char* argv[])
{
    Tube::enable_syntheses();

    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd fullSep(phi,X0);
    double epsilon = 0.1;


    beginDrawing();

    VIBesFigMap fig_map("Example 2");
    fig_map.set_properties(50,50,800,306);



    IntervalVector proj(1);
    proj[0] = Interval(0,8);
    SepProj sepProj(fullSep,proj,epsilon);
    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"g[g]");
    fig_map.axis_limits(x);
    endDrawing();

    return(0);
    

}
