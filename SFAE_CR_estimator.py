import numpy as np
from scipy.optimize import minimize
import pandas as pd
import argparse


def read(file_name, sep=','):
    matrix = np.genfromtxt(file_name,
                           delimiter=sep,
                           dtype=None,
                           names=True,
                           )
    return matrix


def get_parser():
    parser = argparse.ArgumentParser(description='Estimate the weight of cell components')
    parser.add_argument('--input_file', type=str, default=None, help="input expression file")
    parser.add_argument('--num_component', type=int, default='', help="the number of component to estimate")
    return parser.parse_args()


def objectiveFunc_two(i, table):
    a = i[0]
    # ipdb.set_trace()
    print(a)
    formula = 0
    for i in range(len(table)):
        whole, cyto, nuc = table['whole'][i], table['cyto'][i], table['nuc'][i]
        if whole == 0:
            continue
        if (cyto - nuc) * a + nuc <= 0:
            continue
        if whole > max(cyto, nuc) or whole < min(cyto, nuc):
            continue
        formula += \
            (np.log((a * cyto + (1 - a) * nuc) / whole)) ** 2
    return formula,


def objectiveFunc_four(i, table):
    a, b, c = i
    formula = 0
    print(i)
    for index in range(len(table)):
        whole, cyto, nuc, ins, mem = table['total'][index], table['Cytosol'][index], table['Nucleus'][index], \
                                     table['Insoluable'][index], table['Membrane'][index]

        if whole == 0:
            continue
        if whole > max(cyto, nuc, ins, mem) or whole < min(cyto, nuc, ins, mem):
            continue

        combination = a * cyto + b * nuc + c * ins + (1 - a - b - c) * mem
        if combination <= 0:
            continue
        self_add = (np.log(combination / whole)) ** 2
        formula += self_add
    return formula,


def main():
    args = get_parser()
    table_filename = args.input_file
    num_components = args.num_component
    init = [1/num_components] * (num_components - 1)
    init = np.array(init)
    bounds = [(0., 1.)] * (num_components - 1)
    print("Start Optimization\n")
    if num_components == 4:
        table = pd.read_csv(table_filename)
        result = minimize(fun=objectiveFunc_four,
                          x0=init,
                          args=(table,),
                          jac='2-point',
                          method='trust-constr',
                          bounds=bounds)
        print('[Output]')

        print(f"Cytosol weight : {result['x'][0]:.3}\nNucleus weight : {result['x'][1]:.3}\n"
              f"Insoluable weight : {result['x'][2]:.3}\nMembrane weight : {(1 - result['x'][0] - result['x'][1] - result['x'][2]):.3}")

    elif num_components == 2:
        table = pd.read_csv(table_filename)
        print('[Iteration]')

        result = minimize(fun=objectiveFunc_two,
                          x0=init,
                          args=table,
                          jac='2-point',
                          method='trust-constr',
                          bounds=((0., 1.),),
                          )

        print('[Output]')
        print(f"Cytosol weight : {result['x'][0]:.3}\nNucleus weight : {(1 - result['x'][0]):.3}")


if __name__ == '__main__':
    main()
