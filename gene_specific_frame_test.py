import argparse
from collections import OrderedDict
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import math
import itertools
import os
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqUtils import CodonUsage
from Bio import SeqIO
import statsmodels.api as sm
from scipy import stats
import scipy
import pandas as pd
import numpy as np
import plastid
import matplotlib
import matplotlib.pyplot as plt
matplotlib.style.use("ggplot")
fdrc = sm.stats.multipletests


def create_savedir(path_save, prefix):
    try:
        os.mkdir(path_save + "/" + prefix)
    except:
        pass


def repair_array(x):

    def checknan(z):
        if isinstance(z, float):
            return z
        else:
            try:
                return np.array(z)
            except:
                pass

    x = x.apply(lambda z: checknan(z))
    return x


def return_ind(x):
    x = "_".join(x.split("_")[:3])
    return x


def get_mn_mle_fit(counts):
    mle_p = np.sum(counts, axis=0) / np.sum(counts)
    mle_ll = sum([scipy.stats.multinomial.logpmf(
        counts[x], counts[x].sum(), mle_p) for x in range(len(counts))])
    return mle_p, mle_ll


def stack_genes_reps(x):
    stacks = {}
    for f in x.columns:
        s = np.vstack(x[f].dropna().values)
        s = s[~np.any(s < 10, axis=1)]
        if len(s) > 0:
            stacks[f] = s
        else:
            stacks[f] = np.nan

    return pd.Series(stacks)


def transform(xt, yt, zt):
    # Transform function for trianglizing data
    sf = np.sqrt(3)/2
    x = 1 / 2 * (xt + 2 * yt) / (xt + yt + zt)
    y = sf * xt / (xt + yt + zt)

    coords = list(np.array((x, y)) * 100)

    return coords[0], coords[1]


def get_plot_misc():
    # Triangle boundries in triangle plot
    an = np.vstack([[0, 0, 1], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    anchors = transform(an.T[0], an.T[1], an.T[2])
    return anchors


def visualize_trigplot(result_df, cat, show_traj, show_anno, s):
    # Triangle plot

    anchors = get_plot_misc()
    ctr_x, ctr_y = transform(
        result_df["ctr_f0"], result_df["ctr_f1"], result_df["ctr_f2"])
    alt_x, alt_y = transform(
        result_df["alt_f0"], result_df["alt_f1"], result_df["alt_f2"])

    plt.scatter(ctr_x, ctr_y, color="darkgreen", s=s)
    plt.scatter(alt_x, alt_y, color="darkred", s=s)
    plt.plot(anchors[0], anchors[1], color="black", linestyle="--")

    if show_traj == True:
        plt.plot(([ctr_x, alt_x]), ([ctr_y, alt_y]),
                 color="black", linestyle="--", lw=0.15)

    if show_anno == True:
        for index, row in result_df.iterrows():
            plt.annotate(
                s=index[:-2], xy=(transform(row["alt_f0"], row["alt_f1"], row["alt_f2"])), size=5)

    plt.savefig(path_save + "/" + prefix + "/" +
                prefix + cat + ".png", dpi=150)
    plt.close()


def visualize_p_q(result_df, cat, s):
    # Plot a histogram of P-values (red) and Q-values (blue) in the samples

    plt.hist(result_df["pval"], alpha=0.5, bins=s)
    plt.hist(result_df["qval"], alpha=0.5, bins=s)

    plt.savefig(path_save + "/" + prefix + "/" +
                prefix + cat + "_pq_dist.png", dpi=150)
    plt.close()


def run_lr_test(dframe, cat):
    # Likelihood-ratio test

    result_dict = {}

    ctr = dframe.loc["control"]
    alt = dframe.loc["treated"]

    for col in dframe.columns:
        result_dict[col] = {}

        # Minimum coverage
        mincov = np.min(dframe[col].apply(lambda x: np.min(x.sum(1))))
        ctr_scaled_counts = np.round(((dframe[col].loc["control"]
                                       / dframe[col].loc["control"].sum(1)[::, np.newaxis])
                                      * mincov), 0)

        alt_scaled_counts = np.round(((dframe[col].loc["treated"]
                                       / dframe[col].loc["treated"].sum(1)[::, np.newaxis])
                                      * mincov), 0)

        # Find distributions under MLE for control, treated, and combined
        dist_ctr, logL_ctr = get_mn_mle_fit(np.atleast_2d(ctr_scaled_counts))
        dist_alt, logL_alt = get_mn_mle_fit(np.atleast_2d(alt_scaled_counts))
        dist_null, logL_null = get_mn_mle_fit(np.vstack(
            [np.atleast_2d(ctr_scaled_counts), np.atleast_2d(alt_scaled_counts)]))

        # Calculate test statistic
        gstat = 2*((logL_ctr+logL_alt)-logL_null)

        # Calculate p-value from test statistic
        p_val = scipy.stats.chi2.sf(gstat, 2)

        # Calculate effect size Cramer's V - tends to 0 as P-value approaches 1
        effect_size = np.sqrt(
            gstat/(np.sum(np.vstack([np.atleast_2d(ctr[col]), np.atleast_2d(alt[col])])) * 2))

        # Calculate the difference in distributions of control and treated
        diff = (dist_alt-dist_ctr)

        # Multiply effect size with positive and negative distribution differences in each frame
        diff[diff > 0] = diff[diff > 0] / np.sum(diff[diff > 0]) * effect_size
        diff[diff < 0] = - \
            (diff[diff < 0] / np.sum(diff[diff < 0]) * effect_size)

        # Store values for future CSV output
        result_dict[col]["ctr_f0"] = dist_ctr[0]
        result_dict[col]["ctr_f1"] = dist_ctr[1]
        result_dict[col]["ctr_f2"] = dist_ctr[2]

        result_dict[col]["alt_f0"] = dist_alt[0]
        result_dict[col]["alt_f1"] = dist_alt[1]
        result_dict[col]["alt_f2"] = dist_alt[2]

        result_dict[col]["pval"] = p_val
        result_dict[col]["effect_size"] = effect_size

        result_dict[col]["diff_f0"] = diff[0]
        result_dict[col]["diff_f1"] = diff[1]
        result_dict[col]["diff_f2"] = diff[2]

    result_df = pd.DataFrame.from_dict(result_dict, orient="index")
    # If p-value is 1, effect size approaches square root of 0 and is returned NaN. Fill NaN with 0
    result_df.fillna(0, inplace=True)
    result_df["qval"] = fdrc(result_df["pval"], alpha=0.1, method="fdr_bh")[1]
    result_df.to_csv(path_save + "/" + prefix + "/" + prefix + cat + ".csv")
    return result_df

    # Combine the replicates for control and treated


def load_samples_dataframes(path_control, path_treated):

    result_counts_gene = {}
    result_counts_roi = {}

    for root, subd, files in os.walk(path_control):
        for file in files:
            path = root + "/" + file
            counts_df = pd.read_json(path)
            counts_df[counts_df.select_dtypes("object").columns] = counts_df[counts_df.select_dtypes(
                "object").columns].apply(lambda x: repair_array(x))

            result_counts_gene[return_ind(
                file)] = counts_df["Frame_counts"].to_dict()
            result_counts_gene[return_ind(file)]["type"] = "control"

            result_counts_roi[return_ind(file)] = {}
            for key in roi_feats:
                for s in ["_a", "_p", "_e"]:
                    result_counts_roi[return_ind(
                        file)][key+s] = counts_df[key+s].sum()

            result_counts_roi[return_ind(file)]["type"] = "control"

    for root, subd, files in os.walk(path_treated):
        for file in files:
            path = root + "/" + file
            counts_df = pd.read_json(path)
            counts_df[counts_df.select_dtypes("object").columns] = counts_df[counts_df.select_dtypes(
                "object").columns].apply(lambda x: repair_array(x))

            result_counts_gene[return_ind(
                file)] = counts_df["Frame_counts"].to_dict()
            result_counts_gene[return_ind(file)]["type"] = "treated"

            result_counts_roi[return_ind(file)] = {}
            for key in roi_feats:
                for s in ["_a", "_p", "_e"]:
                    result_counts_roi[return_ind(
                        file)][key+s] = counts_df[key+s].sum()
            result_counts_roi[return_ind(file)]["type"] = "treated"

    # Read the dict as dataframe and drop genes without sufficient coverage in all replicates
    df_counts_gene = pd.DataFrame.from_dict(result_counts_gene, orient="index")
    df_counts_roi = pd.DataFrame.from_dict(result_counts_roi, orient="index")

    df_counts_gene = df_counts_gene.dropna("columns").groupby(
        "type").apply(lambda x: stack_genes_reps(x)).dropna("columns")
    df_counts_roi = df_counts_roi.dropna("columns").groupby(
        "type").apply(lambda x: stack_genes_reps(x)).dropna("columns")

    return df_counts_gene, df_counts_roi


def output_results(df_counts_gene, df_counts_roi):
    # Perform the test for each codon and amino acid at A, P and E site; and genes
    for s in ["_a", "_p", "_e"]:
        try:
            features = [x+s for x in codon_dict.keys()]
            result_df = run_lr_test(df_counts_roi[features], "_codon" + s)
            visualize_trigplot(result_df, "_codon" + s, True, False, 5)
        except:
            pass

    for s in ["_a", "_p", "_e"]:
        try:
            features = [x+s for x in aa_dict.keys()]
            result_df = run_lr_test(df_counts_roi[features], "_amino_acid" + s)
            visualize_trigplot(result_df, "_amino_acid" + s, True, True, 5)
        except:
            pass

    result_df = run_lr_test(df_counts_gene, "_gene")
    visualize_trigplot(result_df, "_gene", False, False, 1)
    visualize_p_q(result_df, "_gene", 100)


if __name__ == "__main__":

    aa_dict = dict(CodonUsage.SynonymousCodons)
    codon_dict = dict(CodonUsage.CodonsDict)
    roi_feats = list([*codon_dict.keys(), *aa_dict.keys()])

    parser = argparse.ArgumentParser()
    parser.add_argument('--control')
    parser.add_argument('--treated')
    parser.add_argument('--savepath')
    args = parser.parse_args()
    prefix = args.control.split("/")[-1] + "_" + args.treated.split("/")[-1]
    path_save = args.savepath

    create_savedir(args.savepath, prefix)

    df_counts_gene, df_counts_roi = load_samples_dataframes(
        args.control, args.treated)
    output_results(df_counts_gene, df_counts_roi)
