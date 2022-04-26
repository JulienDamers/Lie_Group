from codac import *
from codac.unsupported import *
from lie_group_ex4_separator import lie_group_ex4_separator

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


def example_4_continuous():
    print("Running example 4\n")
    # Generate reference
    domain = Interval(0, 15)
    timestep = 0.001
    x0 = IntervalVector(4)
    x0.init(Interval(0, 0))
    f = Function("x", "y", "z", "x4", "(cos(z); sin(z); sin(0.4*x4); 1)")
    ctc_lohner = CtcLohner(f)
    a = TubeVector(domain, timestep, f.nb_arg())
    a.set(x0, 0.)
    ctc_lohner.contract(a)
    print("Reference generated ", a)

    epsilon = 0.01

    X0 = IntervalVector(3)
    X0[0] = Interval(-0.1, 0.1)
    X0[1] = Interval(-0.1, 0.1)
    X0[2] = Interval(-0.4, 0.4)

    fullSep = lie_group_ex4_separator(a, X0)

    proj = IntervalVector([[-6, 6], [-15, 0]])
    sep_proj = SepProj(fullSep, proj, epsilon)

    mapping = IntervalVector([[-4, 4], [-4, 4]])

    beginDrawing()
    fig = VIBesFigMap("example_4_continuous_py")
    fig.set_properties(50, 50, 800, 800)
    background_box = IntervalVector(mapping)
    fig.axis_limits(background_box)
    start = time.time()
    SIVIA(mapping, sep_proj, epsilon)
    stop = time.time()
    print("elapsed time test-case 4 continuous: ", stop - start, "s")


if __name__ == "__main__":
    example_4_continuous()

