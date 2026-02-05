import numpy as np
import scipy 
import python_general_util_func as gutil

def main():
    h = 1
    L = 2
    Deq = gutil.computeDeq(h,L)
    print(str(Deq) + f", expected output: 1.5234\n")

    NPR = 3
    Mj = gutil.computeMj(NPR)
    Ma = gutil.computeMa(NPR)
    print(str(Mj) + f", Expected value: 1.3578\n")
    print(str(Ma) + f", Expected value: 1.1606\n")

    f_scr = gutil.calcSt2fND(0.2,Ma,Deq)
    print(rf"f_scr is {f_scr}, expected value: 0.2819")
if __name__ == "__main__":
    main()
