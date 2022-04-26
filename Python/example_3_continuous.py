from codac import *
from codac.unsupported import *

import time


class CtcCn(Ctc):

    def __init__(self, cn, box_to_contract_cn):
        Ctc.__init__(self, len(box_to_contract_cn))

        self.cn = cn
        self.box_to_contract_cn = box_to_contract_cn

    def contract(self, box):
        self.cn.reset_interm_vars()
        self.cn.contract({self.box_to_contract_cn: box})

        return box


def example_3_continuous():
    print("Running example 3\n")
    # Generating the reference
    t_domain = Interval(0, 6)  # Variable à modifier pour changer le temps d'intégration
    time_step = 0.005
    a = TubeVector(t_domain, time_step, 2)
    x0 = IntervalVector([[0.5, 0.5], [0, 0]])  # Initial condition for the reference
    a.set(x0, 0.)
    f = Function("x", "y", "(-x^3-x*y^2+x-y; -y^3-x^2*y+x+y)")
    ctc_lohner = CtcLohner(f)
    ctc_lohner.contract(a)

    print("Reference generated", a)

    epsilon = 0.01  # SIVIA accuracy

    psi = Function("x10", "x20", "t", "r1", "r2",
                   "((x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))-x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)))-r1;\
                   (x20*(0.5*(sqrt(1/(3*exp(-2*t)+1))*cos(t)))+x10*(0.5*(sqrt(1/(3*exp(-2*t)+1))*sin(t))))/(0.25*sqrt(  ((1-(1/(3*exp(-2*t)+1)))/0.75)+( ((1/(3*exp(-2*t)+1))*(4+(1/0.75)))-(4/3) ) * (x10^2+x20^2)  ))-r2)")

    f_dom = Function("x10", "x20", "t", "r1", "r2",
                     "((1-(1/(3*exp(-2*t)+1)))/0.75)+(((1/(3*exp(-2*t)+1))/0.25)-((1-(1/(3*exp(-2*t)+1)))/0.75))*(x10^2+x20^2)")

    # Création des contracteur sur nos fonctions
    ctc_psi = CtcFwdBwd(psi)
    ctc_psi_int = CtcFwdBwd(psi) | CtcFwdBwd(f_dom, Interval.NEG_REALS)

    f_circle = Function("i1", "i2", "( (i1-1.5)^2 + (i2-1.5)^2)")
    ctc_initial = CtcFwdBwd(f_circle, Interval(0, 0.04))
    ctc_not_initial = CtcNotIn(f_circle, Interval(0, 0.04))

    # contracteur extérieur
    cn_out = ContractorNetwork()
    box_out = IntervalVectorVar(3)
    x_init_out = cn_out.create_interm_var(IntervalVector(2))  # x0
    cn_out.add(ctc_psi, [box_out, x_init_out])
    cn_out.add(ctc_initial, [x_init_out])  # x_init in circle
    ctc_cn_out = CtcCn(cn_out, box_out)

    # contracteur intérieur
    cn_in = ContractorNetwork()
    box_in = IntervalVectorVar(3)
    x_init_in = cn_in.create_interm_var(IntervalVector(2))  # x0
    cn_in.add(ctc_psi_int, [box_in, x_init_in])
    cn_in.add(ctc_not_initial, [x_init_in])  # x_init not in circle
    ctc_cn_in = CtcCn(cn_in, box_in)

    sep = SepCtcPair(ctc_cn_in, ctc_cn_out)
    sep_proj = SepProj(sep, -t_domain, 0.01)

    mapping = IntervalVector([[-2, 2], [-2, 2]])

    beginDrawing()
    fig = VIBesFigMap('example_3_continous_py')
    fig.set_properties(50, 50, 800, 800)
    backgound_box = IntervalVector(mapping)
    fig.axis_limits(backgound_box)
    start = time.time()
    SIVIA(mapping, sep_proj, epsilon)
    stop = time.time()
    print("elapsed time test-case 3 continuous: ", stop - start, "s")

if __name__ == "__main__":
    example_3_continuous()

