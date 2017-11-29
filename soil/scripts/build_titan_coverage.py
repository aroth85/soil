from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

import pandas as pd

pandas2ri.activate()


def main(args):
    importr('HMMcopy')
    titan = importr('TitanCNA')

    if args.target_bed_file is None:
        df = titan.correctReadDepth(
            args.tumour_wig_file,
            args.normal_wig_file,
            args.gc_wig_file,
            args.mappability_wig_file,
        )

    else:
        target_df = pd.read_csv(args.target_bed_file, header=None, sep='\t')

        df = titan.correctReadDepth(
            args.tumour_wig_file,
            args.normal_wig_file,
            args.gc_wig_file,
            args.mappability_wig_file,
            targetedSequence=pandas2ri.py2ri(target_df)
        )

    df = pandas2ri.ri2py_dataframe(df)

    df.to_csv(args.out_file, index=False, sep='\t')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--normal-wig-file', required=True)

    parser.add_argument('-t', '--tumour-wig-file', required=True)

    parser.add_argument('-g', '--gc-wig-file', required=True)

    parser.add_argument('-m', '--mappability-wig-file', required=True)

    parser.add_argument('-o', '--out-file', required=True)

    parser.add_argument('--target-bed-file', default=None)

    args = parser.parse_args()

    main(args)
