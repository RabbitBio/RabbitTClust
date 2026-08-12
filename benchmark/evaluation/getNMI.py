#!/usr/bin/env python3
import argparse

import numpy as np
from sklearn.metrics import normalized_mutual_info_score


def load_labels(label_file):
    labels = np.loadtxt(label_file, dtype=np.int64, ndmin=2)
    if labels.shape[0] != 2:
        raise ValueError(
            f"{label_file} must contain exactly two rows: "
            "predicted labels followed by ground-truth labels"
        )
    if labels.shape[1] == 0:
        raise ValueError(f"{label_file} does not contain any labels")
    return labels


def getNMI(arr):
    prediction = arr[0, :]
    ground_truth = arr[1, :]
    return normalized_mutual_info_score(prediction, ground_truth)


def main():
    parser = argparse.ArgumentParser(
        description="Compute NMI from a two-row predicted/ground-truth label file."
    )
    parser.add_argument(
        "label_file",
        help="two-row file containing predicted labels and ground-truth labels",
    )
    args = parser.parse_args()

    labels = load_labels(args.label_file)
    print(f"{args.label_file}:\t{getNMI(labels)}")


if __name__ == "__main__":
    main()
