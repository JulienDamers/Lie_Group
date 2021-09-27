//
// Created by julien-damers on 23/08/2021.
//


#include "codac.h"
#include "codac-rob.h"
#include "codac-capd.h"
#include "tools.h"
#include "chrono"
#include "ctc_cn.h"



using namespace std;
using namespace codac;
using namespace ibex;
using namespace vibes;
using namespace pyibex;

void example_2(bool lohner_done, bool capd_done)
{
    cout << "##########################" << endl;
    cout << "####### Particles ########" << endl;
    cout << "##########################" << endl << endl;

    auto start = chrono::steady_clock::now();
    auto stop = chrono::steady_clock::now();
    codac::CtcEval ctc_eval;

    Interval domain_2_x(0, 12); // Define domain of work on which we want to integrate
    double timestep_2 = 0.1;
    IntervalVector a0_2({{0.0, 0.0},
                         {0.0, 0.0}}); // Define initial condition
    Function f_2("x", "y", "(1;sin(x))"); // Evolution function to integrate

    IntervalVector X0({{0,1},{0,1}});

    TubeVector a_capd_2 = CAPD_integrateODE(domain_2_x, f_2, a0_2, timestep_2);



    IntervalVector frame({{-1,10},{-1,4}});
    ipegenerator::Figure fig(frame,400,153);
    fig.set_graduation_parameters(-1,0.5,-1,0.5);
    fig.set_number_digits_axis_x(1);
    fig.set_number_digits_axis_y(1);


    fig.draw_box(X0,"colorBlindMaybeStroke","colorBlindMaybeStroke");

    double rx,ry;
    IntervalVector particle(2);
    IntervalVector particle_image(2);

    for (int i=0; i<1000; i++)
    {
        rx = ((double) rand() / (RAND_MAX));
        ry = ((double) rand() / (RAND_MAX));
        double x = X0[0].lb()+rx*X0[0].diam();
        double y = X0[0].lb()+ry*X0[1].diam();
        particle[0] = Interval(x);
        particle[1] = Interval(y);

        particle_image[0] = particle[0] + a_capd_2(8.)[0];
        particle_image[1] = particle[1] + a_capd_2(particle_image[0])[1] - a_capd_2(particle[0])[1];

        particle.inflate(0.01);
        particle_image.inflate(0.01);
        fig.draw_box(particle, "colorBlindOutStroke","colorBlindOutFill");
        fig.draw_box(particle_image, "colorBlindInStroke","colorBlindInFill");

    }

    fig.set_color_stroke("black");
    fig.set_color_fill("black");
    fig.set_color_type(ipegenerator::STROKE_AND_FILL);
    fig.draw_tubeVector(&a_capd_2,0,1);

    fig.draw_axis("x1","x2");
    fig.save_ipe("particles.ipe");
    fig.save_pdf("particles.pdf");
}

int main(int argc, char* argv[])
{
    example_2(true,true);
}