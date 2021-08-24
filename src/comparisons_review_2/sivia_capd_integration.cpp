//
// Created by julien-damers on 23/08/2021.
//

//
// Created by julien-damers on 22/07/2021.
//

#include "capd/capdlib.h"
#include "codac.h"
#include "codac-rob.h"
#include "codac-capd.h"
#include "tools.h"
#include "chrono"


using namespace std;
using namespace codac;
using namespace ibex;
using namespace vibes;
using namespace pyibex;
using namespace capd;


IntervalVector integrate(IntervalVector& box)
{
    try{
        // This is vector field for the Rossler system
        IMap vectorField("var:x,y;fun:-1,-sin(x);");
        // set chaotic parameter values
        // the solver, is uses high order enclosure method to verify the existence
        // of the solution. The order is set to 20.
        IOdeSolver solver(vectorField,5);
        solver.setAbsoluteTolerance(1e-10);
        solver.setRelativeTolerance(1e-10);
        ITimeMap timeMap(solver);
        // this is a good approximation of a periodic point
        IVector c(2);
        c[0] = capd::Interval(box[0].lb(),box[0].ub());
        c[1] = capd::Interval(box[1].lb(),box[1].ub());;
        // define a doubleton representation of the interval vector c
        C0HORect2Set s(c);
        // we integrate the set s over the time T
        capd::interval T(8);

        IVector result = timeMap(T,s);
        return ( IntervalVector({{result[0].leftBound(),result[0].rightBound()},
                                 {result[1].leftBound(),result[1].rightBound()}}));
    }catch(exception& e)
    {
        cout << "\n\nException caught!\n" << e.what() << endl << endl;
        return(box);
    }
}




int main(int argc, char* argv[])
{

    cout.precision(12);
    IntervalVector x({{-1,10},{-1,3.2}});
    IntervalVector x0({{0,1},{0,1}});
    double epsilon = 0.1;

    int bisections = 0;

    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,400,153);
    fig.set_graduation_parameters(-1,0.5,-1,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    stack<IntervalVector> s;
    s.push(x);

    auto start = chrono::steady_clock::now();
    while(!s.empty())
    {
        IntervalVector box = s.top();
        s.pop();
        IntervalVector boxOut = integrate(box);

        if (boxOut.is_subset(x0))
        {
            fig.draw_box(box,"colorBlindInStroke","colorBlindInFill");
        }
        else if ( (boxOut & x0)[0].is_empty() || (boxOut & x0)[1].is_empty() )
        {
            fig.draw_box(box,"colorBlindOutStroke","colorBlindOutFill");
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
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    cout << "bissections: " <<bisections <<endl;
    fig.set_opacity(30);
    fig.draw_box(x0,"green","green");
    fig.set_opacity(100);
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("\\Huge\\mathbb{X}_0",0.4,0.5,true);
    fig.set_size_axis_graduation(18.0);
    fig.draw_axis("x1","x2");
    fig.save_ipe("sivia_capd_ex_2.ipe");
    fig.save_pdf("sivia_capd_ex_2.pdf");
}