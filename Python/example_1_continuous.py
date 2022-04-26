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


def example_1_continuous():
    print("Running example 1\n")
    # Generate reference
    domain = Interval(0, 5)
    timestep = 0.01
    x0 = IntervalVector([[0., 0.], [1., 1.]])
    f = Function("x", "y", "(1;-y)");
    ctc_lohner = CtcLohner(f)
    a = TubeVector(domain, timestep, f.nb_arg())
    a.set(x0, 0.)
    ctc_lohner.contract(a)
    print("Reference generated ", a)

    epsilon = timestep

    X0 = IntervalVector([[0, 1], [2, 3]])
    mapping = IntervalVector([[-0.1, 6.5], [-0.2, 3.5]])
    phi = Function("x1", "x2", "t", "(x1+t;x2/exp(t))")
    fullSep = SepFwdBwd(phi, X0)

    sep_proj = SepProj(fullSep, -domain, epsilon)

    beginDrawing()

    fig = VIBesFigMap("Example 1 continuous")
    fig.set_properties(50, 50, 800, 464);
    backgound_box = IntervalVector(mapping)
    fig.axis_limits(backgound_box)
    start = time.time()
    SIVIA(mapping, sep_proj,epsilon)
    stop = time.time()
    print("elapsed time test-case 1 continuous: ", stop - start, "s")


if __name__ == "__main__":
    example_1_continuous()

