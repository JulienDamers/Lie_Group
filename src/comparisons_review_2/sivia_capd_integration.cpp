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


void sivia_capd_integration()
{
    cout.precision(12);
    IntervalVector x({{-1,10},{-1,3.2}});
    IntervalVector x0({{0,1},{0,1}});
    double epsilon = 0.1;

    int bisections = 0;

    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,150,69);
    fig.set_graduation_parameters(-1,1,-1,1);
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
            fig.set_current_layer("inner");
            fig.draw_box(box,"grey","white");
        }
        else if ( (boxOut & x0)[0].is_empty() || (boxOut & x0)[1].is_empty() )
        {
            fig.set_current_layer("outer");
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
                fig.set_current_layer("uncertain");
                fig.draw_box(box,"colorBlindInStroke","colorBlindInFill");

            }
        }

    }
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    cout << "bissections: " <<bisections <<endl;
    fig.set_opacity(75);
    fig.draw_box(x0,"black","white");
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

void comparison_sivia_capd()
{
    IntervalVector X0({{0,1},{0,1}});
    IntervalVector x({{-1,10},{-1,3.2}});
    ibex::Function phi("x1","x2","t","(x1-t;x2+cos(x1)-cos(x1-t) )");
    SepFwdBwd* fullSep;
    fullSep = new SepFwdBwd(phi,X0);
    double epsilon = 0.1;

    vector<float> projection_times{8.};

    vector<Sep*> seps;
    vector<IntervalVector*> projections;
    for (size_t i=0;i<projection_times.size();i++)
    {
        IntervalVector *proj = new IntervalVector(1); // Defining interval to project
        (*proj)[0] = Interval(projection_times[i],projection_times[i]);
        projections.push_back(proj);
        SepProj *sepProj = new SepProj(*fullSep,*proj,epsilon);
        seps.push_back(sepProj);

    }
    Array<Sep> ar_sep(seps);
    SepUnion usep (ar_sep);

    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,150,69);
    fig.set_graduation_parameters(-1,1,-1,1);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    auto start = chrono::steady_clock::now();
    sivia_article(x,usep,epsilon,fig);
    auto stop = chrono::steady_clock::now();
    cout << "elapsed time: " << chrono::duration_cast<chrono::milliseconds>(stop - start).count() << " ms" <<endl;
    fig.set_opacity(75);
    fig.draw_box(X0,"black","white");
    fig.set_opacity(100);
    fig.add_layer("text");
    fig.set_current_layer("text");
    fig.set_color_stroke("black");
    fig.set_color_type(ipegenerator::STROKE_ONLY);
    fig.draw_text("{\\Huge$\\mathbb{X}_0$}",0.45,0.45,false);
    fig.draw_text("{\\Huge$\\mathbb{X}_8$}",8.45,1.9,false);
    fig.set_size_axis_graduation(18.0);
    fig.draw_axis("x1","x2");
    fig.save_ipe("comparison_sivia_capd.ipe");
    fig.save_pdf("comparison_sivia_capd.pdf");

    int n = seps.size();
    for (int i =0; i<n; i++)
    {
        delete(projections[i]);
        delete(seps[i]);
    }
    return;
}


int main(int argc, char* argv[])
{
    sivia_capd_integration();
    Tube::enable_syntheses();
    comparison_sivia_capd();

}