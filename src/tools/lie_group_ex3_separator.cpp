//
// Created by julien-damers on 18/12/2019.
//

#include "lie_group_ex3_separator.h"

using namespace codac;
using namespace ibex;
using namespace std;


lie_group_ex3_separator::lie_group_ex3_separator(TubeVector* a, IntervalVector* constraint) : Sep(3)
{

    reference = a;
    constraintVector = constraint;
    constraintDomain = Interval(0,POS_INFINITY);


    fy = new ibex::Function("x1","x2","z1","z2",
                            "((sqrt(3)*(x1*z1-x2*z2))/sqrt(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2)); \
                                 (sqrt(3)*(x2*z1+x1*z2))/sqrt(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2)))");
    fdom = new ibex::Function("x1","x2","z1","z2","(1-(z1^2+z2^2)+(4*(z1^2+z2^2)-1)*(x1^2+x2^2))");

    sepPhi = new ibex::SepFwdBwd(*fy,*constraintVector);
    sepDom = new ibex::SepFwdBwd(*fdom,constraintDomain);
    fullSep = new ibex::SepInter(*sepDom,*sepPhi);

}


void lie_group_ex3_separator::separate(ibex::IntervalVector& Xin,ibex::IntervalVector& Xout )
{

    Interval t = Xin[2];

    IntervalVector Xinit(Xin);
    IntervalVector z(2);
    z.init(Interval(NEG_INFINITY,POS_INFINITY));

    if (t.is_empty())
    {
        Xout[0].set_empty();
        Xout[1].set_empty();
        Xout[2].set_empty();
        cout << "z empty" << endl;
        return;
    }


    z = z & (*reference)(t);
    if (z[0].is_empty()||z[1].is_empty())
    {
        Xout[0].set_empty();
        Xout[1].set_empty();
        Xout[2].set_empty();
        cout << "z empty" << endl;
        return;
    }


    IntervalVector fullBoxIn(4);
    fullBoxIn[0] = Xin[0];
    fullBoxIn[1] = Xin[1];
    fullBoxIn[2] = z[0];
    fullBoxIn[3] = z[1];



    IntervalVector fullBoxOut(fullBoxIn);
    IntervalVector fullBoxSave(fullBoxIn);
    IntervalVector zout(z);
    Interval tout(t);

    fullSep->separate(fullBoxIn,fullBoxOut);

    Xin[0] = fullBoxIn[0];
    Xin[1] = fullBoxIn[1];
    z[0] = fullBoxIn[2];
    z[1] = fullBoxIn[3];

    Xout[0] = fullBoxOut[0];
    Xout[1] = fullBoxOut[1];
    zout[0] = fullBoxOut[2];
    zout[1] = fullBoxOut[3];


    t = t & reference->invert(z,t);
    tout = tout & reference->invert(zout,tout);


    assert((t|tout)==Xinit[2]);

    if(t.is_empty())
    {
        Xin[0].set_empty();
        Xin[1].set_empty();
        Xin[2].set_empty();
        Xout[2] = tout;
        return;

    }


    if (tout.is_empty())
    {
        Xout[0].set_empty();
        Xout[1].set_empty();
        Xout[2].set_empty();
        Xin[2] = t;
        return;
    }


    Xin[2] = t;
    Xout[2] = tout;


}


lie_group_ex3_separator::~lie_group_ex3_separator()
{

    delete sepPhi;
    delete fy;
    delete fdom;
    delete sepDom;
    delete fullSep;

}


