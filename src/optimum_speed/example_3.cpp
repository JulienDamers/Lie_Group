//
// Created by julien-damers on 16/12/2019.
//

#include <iostream>
#include <cstdlib>
#include "tools.h"
#include "tubex.h"
#include "tubex-capd.h"
#include "chrono"



using namespace std;
using namespace tubex;
using namespace ibex;
using namespace vibes;
using namespace pyibex;





int main(int argc, char* argv[])
{
    Tube::enable_syntheses();


    Interval domain(0,4);
    double timestep = 0.01;
    IntervalVector x0(2);
    x0[0] = Interval(0.5,0.5);
    x0[1] = Interval(0,0);
    // Generate the reference  as we do not have a analytical expression for it
    Function f("x","y","(x^3+x*y^2-x+y; y^3+x^2*y-x-y)");
    // CAPD integration version
    TubeVector a = CAPD_integrateODE(domain,f,x0,timestep);
    cout << "Reference generated " << a << endl;

    double epsilon = 0.01; // precision of the sivia

    IntervalVector X0({{0.4,0.6},{-0.1,0.1}}); // Starting box constraint
    IntervalVector x({{-1.2,1.2},{-1.2,1.2}}); // Space to be explored
    lie_group_ex3_separator fullSep(&a,&X0);
    IntervalVector proj(1);
    proj[0] = Interval(0,4);
    SepProj sepProj(fullSep,proj,0.01);





    beginDrawing();
    VIBesFigMap fig("Example 3");
    fig.set_properties(50,50,800,800);

    auto start = chrono::steady_clock::now();
    sivia(x,sepProj,epsilon);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;

    drawCircle(0,0,1);
    drawBox(X0[0].lb(),X0[0].ub(),X0[1].lb(),X0[1].ub(),"#00FF00A6[#00FF00A6]");

    fig.set_tube_max_disp_slices(20000);
    fig.add_tube(&a,"reference",0,1);
    fig.set_tube_color(&a, "k[k]");
    fig.show();
    fig.axis_limits(x);

    endDrawing();
}